/********************************************************************************/
/*                                                                              */
/*  grow.c                                                                      */
/*  functions used to predict values for the CONIFERS growth model              */
/*                                                                              */
/********************************************************************************/

/* 	$Id: grow.c 878 2012-04-04 19:27:41Z hamannj $	 */

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"

static int compare_plants_by_plot_plant(
				 const void *ptr1,
				 const void *ptr2 );

/* local function declarations */
void __stdcall project_plant_list( 
      unsigned long           *return_code,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_points,
      struct PLOT_RECORD      *plots_ptr,
      double		          *x0,
      unsigned long           endemic_mortality,      
      unsigned long           sdi_mortality,          
      int                     hcb_growth_on,          
      unsigned long           use_precip_in_hg,       
      unsigned long           use_rand_err,
      unsigned long		      variant,
   unsigned long            use_genetic_gains,
   unsigned long			plantation_age,
   unsigned long            yrst,
   unsigned long            *n_years_after_planting );


static void project_plot( 
   unsigned long           *return_code,
   unsigned long           n_plants,
   struct PLANT_RECORD     *plants_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   unsigned long           endemic_mortality,      
   int                     hcb_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err,
   unsigned long	        variant,
   struct SUMMARY_RECORD   *before_sums,
   unsigned long           use_genetic_gains,
   unsigned long			plantation_age,
   unsigned long            yrst,
   unsigned long            *n_years_after_planting );


/********************************************************************************/
/* project the plant list                                                       */
/* this function will project each plot for one year                            */
/* to project the entire sample for more than one year, this function           */
/* needs to be called once for each year                                        */
/********************************************************************************/
void __stdcall project_plant_list( 
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   unsigned long           n_plants,
   struct PLANT_RECORD     *plants_ptr,
   unsigned long           n_points,
   struct PLOT_RECORD      *plots_ptr,
   double		           *x0, 		   
   unsigned long           endemic_mortality,      
   unsigned long           sdi_mortality,          
   int                     hcb_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err,
   unsigned long	        variant,
   unsigned long            use_genetic_gains,
   unsigned long			plantation_age,
   unsigned long            yrst,
   unsigned long            *n_years_after_planting )
{
   unsigned long           i;
   struct  PLOT_RECORD     *plot_ptr;
   double                  max_sdi;
   struct SUMMARY_RECORD   before_sums;
   struct SUMMARY_RECORD   after_sums;
   double		           mortality_proportion;

   double  plants_removed;
   double  ba_removed;

    qsort(   plants_ptr, 
       n_plants, 
       sizeof( struct PLANT_RECORD ), 
       compare_plants_by_plot_plant ); 


   /* update the plot variabels before you procede */
   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );
   calc_plot_stats_2(   return_code, 
                           n_species,
                           species_ptr,
                           n_coeffs,
                           coeffs_ptr,
                           n_plants,
                           plants_ptr,
                           n_points,
                           plots_ptr );

   if( *return_code != CONIFERS_SUCCESS )
   {
      return;
   }

   /* update the total summaries before the plant list is projected */
   update_total_summaries(  return_code,
                               n_points,
                               n_plants,
                               n_species,
                               species_ptr,
                               n_coeffs,
                               coeffs_ptr,
                               plants_ptr,
                               &before_sums );

   if( *return_code != CONIFERS_SUCCESS )
   {
      return;
   }


   /* for each plot, project it forward one year */
   plot_ptr = &plots_ptr[0];
   for( i = 0; i < n_points; i++, plot_ptr++ )
   {

      project_plot( return_code,
                      n_plants,
                      plants_ptr,
                      plot_ptr,
                      n_species,
                      species_ptr,
                      n_coeffs,
                      coeffs_ptr,
                      endemic_mortality,
                      hcb_growth_on,          
                      use_precip_in_hg,  
                      use_rand_err,
                      variant,    
                      &before_sums,
                      use_genetic_gains,
	                  plantation_age,
	                  yrst,
                      n_years_after_planting );

      if( *return_code != CONIFERS_SUCCESS )
      {
            /* todo: should set some warning in here */
            //* warning = WARNING_TYPE;
	 return;
	 //continue;
      }


   }

   /* now, calculate the limiting sdi mortality value for  */
   /* all the trees in the plant list                      */
    
   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );

   /* calcuate rd */
   /* update the total summaries before the plant list is projected */
   update_total_summaries(  return_code,
                               n_points,
                               n_plants,
                               n_species,
                               species_ptr,
                               n_coeffs,
                               coeffs_ptr,
                               plants_ptr,
                               &after_sums );

   /* need_error_trap_here */
   /* if max sdi limit switch is on, then      */
   /* MOD012 */
   //Rprintf( "if( sdi_mortality ) check == 1 @ %s, %d, then ()\n", __FILE__, __LINE__ );
   if( sdi_mortality )
   {
      calc_max_sdi( return_code,
                      n_species,
                      species_ptr,
                      n_coeffs,
                      coeffs_ptr,
                      n_plants,
                      plants_ptr,
                      n_points,
                      &max_sdi  );

      /* it is choking here. */
      if( *return_code != CONIFERS_SUCCESS )
      {
	 //Rprintf( "calc_max_sdi *return_code != CONIFERS_SUCCESS @ %s, %d\n", __FILE__, __LINE__ );
	 return;
      }


      /* this might have been changed since the code was ported */
      /* to the R interface, but a check needs to happen to	*/
      /* insure that x0 isn't not reset. A check will be placed */
      /* here to verify x0 has not been set.			*/      
      if( *x0 == 0.0 ) 
      {
	 //Rprintf( "%s, %d\n", __FILE__, __LINE__ );
  	 //Rprintf( "*x0 == 0.0 -- calling calc_hann_wang_x0\n" );
  	 //Rprintf( "after_sums.sdi %f\n", after_sums.sdi );
  	 //Rprintf( "max_sdi = %lf\n", max_sdi );

	 /* this seems to be broken */
	 calc_hann_wang_x0( return_code,
	                    after_sums.sdi,
	                    after_sums.bh_expf,
	                    max_sdi,
	                    before_sums.rel_density,
	                    after_sums.rel_density,
	                    x0);
	 
	 if( *return_code != CONIFERS_SUCCESS )
	 {
	    //Rprintf( "*return_code != CONIFERS_SUCCESS, %s, %d\n", __FILE__, __LINE__ );
	    //Rprintf( "*x0 == 0.0 -- calling calc_hann_wang_x0\n" );
	    return;
	 }
      }
	
      if( *x0 > 0.0 )
      {
	 //Rprintf( "*x0 > 0.0 -- calling calc_sdi_mortality\n" );

	 calc_sdi_mortality( return_code,
	                    after_sums.qmd,
	                    after_sums.sdi,
	                    max_sdi,
	                    *x0,
	                    after_sums.bh_expf,
	                    &mortality_proportion);
	 if( *return_code != CONIFERS_SUCCESS )
	 {
	    return;
	 }

	 /*    this next section does the sdi mortality   */
	 plot_ptr = &plots_ptr[0];
	 for( i = 0; i < n_points; i++, plot_ptr++ )
	 {
	    thin_plot(  return_code,
	                n_plants,
	                plants_ptr,
	                plot_ptr,
	                n_species,
	                species_ptr,
	                n_coeffs,
	                coeffs_ptr,
	                //NULL,
	                0,
	                DO_SDI_MORT,
	                mortality_proportion, 
	                &plants_removed, 
	                &ba_removed );
            
	    if( *return_code != CONIFERS_SUCCESS )
	    {
	       continue;/* error trap here */
	    }
	 }
      }
   } /* end of the sdi_mortality switch */


    (*n_years_after_planting)++;

}








/********************************************************************************/
/* project_plot                                                                 */
/********************************************************************************/
/*  Description :   projects all the plants on the plot forward one year        */
/*                  into the future                                             */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   September 12, 1999                                          */
/*  Returns     :   void                                                        */
/*  Comments    :   this function updates the following variables for each		*/
/* plant on the plot that is projected											*/
/*      plant_ptr->d6           += plant_ptr->d6_growth;						*/
/*      plant_ptr->d12          += plant_ptr->d12_growth;						*/
/*      plant_ptr->dbh          += plant_ptr->dbh_growth;						*/
/*      plant_ptr->tht          += plant_ptr->tht_growth;						*/
/*      plant_ptr->cr           += plant_ptr->cr_growth;						*/
/*      plant_ptr->crown_width  += plant_ptr->cw_growth;						*/
/*      plant_ptr->expf         -= plant_ptr->expf_change;						*/
/*  Arguments   :   unsigned long *return_code  - pointer to a return code      */
/*                  unsigned long n_plants      - total number fo plants in the */
/*                      plants pointer array                                    */
/*                  struct PLANT_RECORD *plants_ptr - array of plants in the    */
/*                      sample to be projected                                  */
/*                  struct PLOT_RECORD  *plot_ptr   - pointer to the current    */
/*                      plot that is to be grown                                */
/*                  unsigned long n_species - size of the species_ptr           */
/*                  struct SPECIES_RECORD   *species_ptr - array of             */
/*                      SPECIES_RECORD's that hold species specific information */
/*                  unsigned long   n_coeffs - sizes of the coeffs_ptr array    */
/*                  struct COEFFS_RECORD *coeffs_ptr - array of coefficients    */
/*                      that are used to project the individual plants on the   */
/*                      plot.                                                   */
/********************************************************************************/
static void project_plot( 
   unsigned long           *return_code,
   unsigned long           n_plants,
   struct PLANT_RECORD     *plants_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   unsigned long           endemic_mortality,      
   int                     hcb_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err,
   unsigned long	       variant,
   struct SUMMARY_RECORD   *before_sums,
   unsigned long           use_genetic_gains,
   unsigned long			plantation_age,
   unsigned long            yrst,
   unsigned long            *n_years_projected )
{

   unsigned long   i;
   struct  PLANT_RECORD    *plant_ptr;
   double                  si_30;
   double                  h40;
   unsigned long           genetics_age_cut;


   genetics_age_cut=0;
   /* for each plot, project it forward one year */
   if(variant == CONIFERS_SMC)
   {
      h40   =  before_sums->height_40;
      si_30 =  plot_ptr->site_30;
      get_age_cut(return_code,
                  h40,
                  si_30,
                  &genetics_age_cut);
   }


    /* loop over the trees and if the plant is on the current plot  */
    /* then project the plant forward for one year                   */
   plant_ptr = &plants_ptr[0];
   for( i = 0; i < n_plants; i++, plant_ptr++ )
   {

      /* if you're not on the plot, don't grow it */
      /* you should really be using the function  */
      /* that returns the index values for the    */
      /* trees on the plot to get the thing done  */
      /* just iterate through the plant list and  */
      /* pitch the plant if it's not on the plot  */
      if( plant_ptr->plot != plot_ptr->plot )
      {
	    continue;
      }

      /* for each tree on the plot....    */
      switch( variant )
      {
	    case CONIFERS_SWO:
	      swo_project_plant(return_code,
			                n_species, 
			                species_ptr, 
			                n_coeffs, 
			                coeffs_ptr,
			                plant_ptr, 
			                plot_ptr,
			                endemic_mortality,  
			                hcb_growth_on,      
			                use_precip_in_hg,    
			                use_rand_err  );  
		break;

	    case CONIFERS_SMC:
	      smc_project_plant(return_code,
			                n_species, 
			                species_ptr, 
			                n_coeffs, 
			                coeffs_ptr,
			                plant_ptr, 
			                plot_ptr,
			                endemic_mortality,  
			                hcb_growth_on,      
			                use_precip_in_hg,    
			                use_rand_err,
			                before_sums,
                            use_genetic_gains,
                            genetics_age_cut);   
		break;

		/* todo: add the code to project each plant for your variant */
	    case CONIFERS_SWOHYBRID:
			swo_hybrid_project_plant( return_code,
			                n_species, 
			                species_ptr, 
			                n_coeffs, 
			                coeffs_ptr,
			                plant_ptr, 
			                plot_ptr,
			                endemic_mortality,  
			                hcb_growth_on,      
			                use_precip_in_hg,    
			                use_rand_err );
			break;

		
		/* todo: add the code to project each plant for your variant */
	    case CONIFERS_CIPS:
	      cips_project_plant(return_code,
			                n_species, 
			                species_ptr, 
			                n_coeffs, 
			                coeffs_ptr,
			                plant_ptr, 
			                plot_ptr,
			                endemic_mortality,  
			                hcb_growth_on,      
			                use_precip_in_hg,    
			                use_rand_err,
			                before_sums,
                            use_genetic_gains,
                            genetics_age_cut,
							plantation_age,
							yrst,
                            n_years_projected ); 
          

        /* todo: put a loop in here and update the dbh using the */
          /* check in project_plant to impute the missing dbh   */
          /* values and make sure you get the end of cycle veg cover    */


		break;


	    default:
	      swo_project_plant(return_code,
			                n_species, 
			                species_ptr, 
			                n_coeffs, 
			                coeffs_ptr,
			                plant_ptr, 
			                plot_ptr,
			                endemic_mortality,  
			                hcb_growth_on,      
			                use_precip_in_hg,    
			                use_rand_err  );     
		break;
      }

      if( *return_code != CONIFERS_SUCCESS )
      {
	     *return_code = FAILED_PROJECT_PLANT;  
	     return;
      }

      
      /* update the current tree values...                      */
      /* todo: this should be moved to a increment/add function */
      plant_ptr->d6           += plant_ptr->d6_growth;
      plant_ptr->d12          += plant_ptr->d12_growth;
      plant_ptr->dbh          += plant_ptr->dbh_growth;
      plant_ptr->tht          += plant_ptr->tht_growth;
      plant_ptr->cr           += plant_ptr->cr_growth;
      plant_ptr->crown_width  += plant_ptr->cw_growth;
      plant_ptr->expf         -= plant_ptr->expf_change;

      /* update the extra data for the tree                     */
      /* TODO: we should probably put this into another function */
      plant_ptr->basal_area    = plant_ptr->dbh * plant_ptr->dbh * FC_I;
      plant_ptr->d6_area       = plant_ptr->d6 * plant_ptr->d6 * FC_I;
      plant_ptr->d12_area       = plant_ptr->d12 * plant_ptr->d12 * FC_I;
      plant_ptr->crown_area    = plant_ptr->crown_width * 
                                    plant_ptr->crown_width * MY_PI / 4.0;
      plant_ptr->pct_cover     = 100.0 * plant_ptr->expf * 
                                    plant_ptr->crown_area/SQ_FT_PER_ACRE;
   
    }





}







static int compare_plants_by_plot_plant(
				 const void *ptr1,
				 const void *ptr2 )
{
  struct PLANT_RECORD   *pt1_ptr;
  struct PLANT_RECORD   *pt2_ptr;

  pt1_ptr = (struct PLANT_RECORD*)ptr1;
  pt2_ptr = (struct PLANT_RECORD*)ptr2;

  if( pt1_ptr->plot < pt2_ptr->plot )
    {
      return -1;
    }
  if( pt1_ptr->plot > pt2_ptr->plot )
    {
      return 1;
    }
  else
    {

      if( pt1_ptr->plant < pt2_ptr->plant )
        {
          return -1;
        }
      if( pt1_ptr->plant > pt2_ptr->plant )
        {
          return 1;
        }
        else
        {
          return 0;
        }

  }


}

