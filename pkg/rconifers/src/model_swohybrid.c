/********************************************************************************/
/*                                                                              */
/*  model.c                                                                     */
/*  functions used to predict values for the CONIFERS growth model              */
/*                                                                              */
/********************************************************************************/
                                                                               
/*------------------------------------------------------------------------------*/
/*                           INDEX OF FUNCTIONS                                 */
/*------------------------------------------------------------------------------*/
/*                                                                              */
/*  Static Functions                     Coefs   Comments                       */
/*                                                                              */
/*  S1 :  calc_crown_width               CW  new new model refit nov09, 1999    */
/*  S2 :  calc_max_crown_width           MW      need to check coeffs nov3      */
/*  S3 :  calc_crown_ratio               CR      refit May   08                 */
/*  S4 :  calc_d6_from_total_height      MS      Correct Oct 99                 */
/*  S5 :  calc_d6_from_ht_and_dbh        DH      Correct Nov 99                 */
/*  S6 :  calc_dbh_from_height           OD      Correct but not needed Nov 99  */
/*  S7 :  calc_dbh_from_height_and_d6    DH      Correct Nov 99                 */
/*  S8 :  calc_exp_from_cover_and_ca             Correct Nov 99                 */
/*  S9 :  calc_height_from_d6            MS      Correct Nov 99                 */
/*  S10:  calc_height_from_dbh           OD      Correct Nov 99                 */
/*  S11:  calc_nstems_from_cw_and_height ST      Correct Jan 00                 */
/*  S12:  calc_cfvol4                    V4      added   Jan 00                 */                 
/*                                                                              */
/*                                                                              */
/*  Dynamic Functions                                                           */
/*                                                                              */
/*  D1:   calc_d6_growth                         Line  786     modified 3/04    */
/*  D2:   calc_dbh_growth                        Line  913     checked          */
/*  D3    calc_height_growth                     Line  986     checked          */
/*  D4    calc_endemic_mort                      Line  XXX     not checked      */
/*  D5    calc_hcb_growth                        Line  XXX     checked          */
/*  D6    calc_cw_growth                                                        */
/*                                                                              */
/*------------------------------------------------------------------------------*/

/* 	$Id: model_swohybrid.c 840 2011-11-22 23:47:23Z hamannj $	 */

/* #ifndef lint */
/* static char vcid[] = "$Id: model_swohybrid.c 840 2011-11-22 23:47:23Z hamannj $"; */
/* #endif /\* lint *\/ */


//#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"

void swo_hybrid_calc_crown_width( 
    unsigned long   *return_code,
    double          d6_area,
    double          total_height,
    double          *pred_crown_width,
    double          *pred_crown_area,
    double          *coeffs_ptr);

void swo_hybrid_calc_max_crown_width( 
    unsigned long   *return_code,
    double          dbh,
    double          total_height,
    double          *pred_max_crown_width,
    double          *coeffs_ptr    );

void swo_hybrid_calc_crown_ratio(
    unsigned long *return_code,  
    double  total_height,
    double  d6,
    double  *pred_cr,
    double  *coeffs_ptr    );

void swo_hybrid_calc_d6_from_total_height(       
    unsigned long   *return_code,
    double          total_height, 
    double          *pred_d6,
    double          *coeffs_ptr    );

void swo_hybrid_calc_d6_from_ht_and_dbh(       
    unsigned long   *return_code,
    double          total_height,
    double          dbh,
    double          *pred_d6,
    double          *coeffs_ptr    );

void swo_hybrid_calc_dbh_from_height(       
    unsigned long   *return_code,
    double          total_height,
    double          *pred_dbh,
    double          *coeffs_ptr    );

void swo_hybrid_calc_dbh_from_height_and_d6(       
    unsigned long   *return_code,
    double          d6,
    double          total_height,
    double          *pred_dbh,
    double          *coeffs_ptr );

void swo_hybrid_calc_exp_from_cover_and_ca(
    unsigned long   *return_code,
    double          pct_cover,
    double          crown_area,
    double          *pred_expf);

void swo_hybrid_calc_height_from_d6(       
    unsigned long   *return_code,
    double          d6,
    double          *pred_total_height,
    double          *coeffs_ptr    );


void swo_hybrid_calc_height_from_dbh(       
    unsigned long   *return_code,
    double          dbh,
    double          *pred_total_height,
    double          *coeffs_ptr    );

void swo_hybrid_calc_nstems_from_cw_and_height(
    unsigned long   *return_code,
    double          total_height,
    double          crown_width,
    long            *pred_n_stems,
    double          *coeffs_ptr      );

void swo_hybrid_calc_volume(
    unsigned long   *return_code,
    double          total_height,
    double          dbh,
    double          *pred_volume,
    double          *coeffs_ptr );

void swo_hybrid_calc_biomass(
    unsigned long   *return_code,
    double          total_height,
    double          basal_diameter,
	double          crown_width,
	double          dbh,
    double     		*pred_biomass,
    double          *coeffs_ptr );


void swo_hybrid_calc_d6_growth(
	unsigned long   *return_code,
    double          height_growth,
    double          crown_width,
    double          total_height,
    double          d6,
    double          cat_shrubs,
    double          cat_conifers,
    double          cat_hardwoods, 
    double          ca_shrubs,
    double          ca_conifers,
    double          ca_hardwoods,
	double          h20_holding_capacity,
    double          *pred_d6_growth,
    double          *coeffs_ptr,
    unsigned long   plant_type);

//	     dghat  = exp(b0
//             + b1 * log(hghat1)
//             + b2 * cw^0.5
//             + b3 * d11
//             + b4 * (catcon)^2
//             + b5 * (cathwsh)^2
//             + b6 * log(SRT))

void swo_hybrid_calc_dbh_growth(   
    unsigned long   *return_code,
    double          total_height,
    double          height_growth,
    double          crown_ratio,
    double          current_dbh,
    double          current_d6,
    double          d6_growth,
    double          *pred_dbh_growth,
    double          *coeffs_ptr );

void swo_hybrid_calc_height_growth(
	unsigned long   *return_code,	
    double          total_height, 
    double          crown_ratio,
    double          d6_area,
    double          precip, 
    double          h20_holding_capacity,
    double          slope,
    double          aspect,
    double          cacon,
    double          catcon,
    double          cahw,
    double          cash,
    double          cathw,
	double          catsh,
	double          basal_d,
    double          elevation,
    double          random_norm_0_1,
    double          random_unif_0_1a,
    double          random_unif_0_1b,
    long            ind_random,
    double          prob_browse,
    double          prob_top_damage,
    unsigned long   use_precip, 
    
	double			min_temp,		/* species specific minumum temperature, in C	*/
	double			max_temp,		/* species specific maximum temperature, in C	*/
	double			opt_temp,		/* species specific optimal temperature, in C	*/
	double			*tday_c,		/* tday_c is a vector of monthly temps, in C	*/
	double			*srad,			/* solar radiation, vector[12]					*/

    double          *height_growth,
    double          *coeffs_ptr,
	unsigned long   plant_type);

void swo_hybrid_calc_cr_growth(
    unsigned long   *return_code,
    int             hcb_growth_on,
    double          total_height,
    double          height_growth,
    double          crown_ratio,
    double          conifer_ca,
    double          hardwood_ca,
    double          shrub_ca,
    double          uniform_0_1,
    double          *cr_growth,
    double          *coeffs_ptr );


void swo_hybrid_calc_cw_growth(   
	unsigned long   *return_code,
	double          total_height,
	double		    height_growth,
	double		    crown_width,
	double		    ca_conifers,
	double		    ca_hardwoods,
	double		    ca_shrubs,
	double		    catcon,
	double		    unif_0_1,
	double          *pred_cw_growth,
	double          *coeffs_ptr,
    unsigned long   plant_type);



void swo_hybrid_calc_endemic_mortality(   
    unsigned long   *return_code,
    double          expansion_factor,
    double          *pred_mortality,
    double          *coeffs_ptr );




/********************************************************************************/
/* swohybrid_impute                                                             */
/********************************************************************************/
/*  Description :   This function fills in the missing values for the plant     */
/*                  list. This function makes two passes. The first pass fills  */
/*                  in the missing dbh,d6, and height then calculated the       */
/*                  plot level stats and plant values in taller variables       */
/*                  before making the second pass to fill in crown ratio.       */
/*                  This is a complete rewrite of fill_in_missing_values.       */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   September 12, 1999                                          */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
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
/*                  unsigned long   variant                                     */
/********************************************************************************/
void swohybrid_impute( 
			    unsigned long           *return_code,
			    unsigned long           n_species,
			    struct SPECIES_RECORD   *species_ptr,
			    unsigned long           n_coeffs,
			    struct COEFFS_RECORD    *coeffs_ptr,
			    unsigned long           variant,
			    unsigned long           n_plants,
			    struct PLANT_RECORD     *plants_ptr,
			    unsigned long           n_points,
			    struct PLOT_RECORD      *plots_ptr,
			    double                  fixed_plot_radius,
			    double                  min_dbh,
			    double                  baf )
{

  unsigned long   i;
  struct  PLANT_RECORD    *plant_ptr = NULL;
  struct  PLOT_RECORD     *plot_ptr = NULL;
  struct  COEFFS_RECORD   *c_ptr = NULL;
  unsigned long   error_count = 0;

//  double  bait[PLANT_TYPES];
//  double  cait[PLANT_TYPES];

//  double        cait_c;
//  double        cait_h;
//  double        cait_s;

  error_count  = 0;
  *return_code = CONIFERS_SUCCESS;

  /* first check to make sure there are plants in the array */
  if( n_plants <= 0 || plants_ptr == NULL ) 
    {
      *return_code = INVALID_PLANT_COUNT;
      return;
    }

  /* first check to make sure there are plants in the array */
  if( n_points <= 0 || plots_ptr == NULL ) 
    {
      *return_code = INVALID_PLOT_COUNT;
      return;
    }

  /* FIRST PASS */
  /* make a first pass to calculate the basic data for the plants */
  /* fill in missing dbh, d6, and total heights */
  plant_ptr = &plants_ptr[0];
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {

      plant_ptr->errors = E_OKDOKEY; /* default value for error is set=ok */ 

      c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

      /* don't go any further if you don't have a functional species code */
      if( c_ptr == NULL )
	{
	  plant_ptr->errors |= E_INVALID_SPECIES;
	  error_count += 1;
	  continue;
	}

      /* only fill in missing values for stocked plots */
      if( is_non_stocked( c_ptr ) )
	{
	  continue;
	}

      /* if the stem is all below d6  */
      if( plant_ptr->tht < 0.50 && !is_non_stocked( c_ptr ) )
	{        
	  plant_ptr->errors |= E_INVALID_HEIGHT;
	  error_count += 1;
	}

      /* if the total height <= 4.5 and there's a dbh obs */
      /* this error triggers on plants that are exactly 4.5 feet tall */
	/* and have a positive dbh observation                          */
	if( plant_ptr->tht < 4.5 && plant_ptr->dbh > 0.0 )
	  {
	    plant_ptr->errors |= E_INVALID_DBH;
	    error_count += 1;
	  }

	/* if it's a tree with a dbh obs */
	if(!is_tree( c_ptr) && plant_ptr->dbh > 0.0)
	  {
	    plant_ptr->errors |= E_INVALID_DBH;
	    error_count += 1;
	  }


	/* if the plant doesn't have a dbh, fill that in    */ 
	/* if the stem is between 6" and 4.5 feet tall      */
	if( plant_ptr->tht >= 0.50 )
	  {
	    if( is_tree( c_ptr ) )
	      {
		if( plant_ptr->d6 == 0.0 )
		  {
		    if(plant_ptr->dbh > 0.0 && plant_ptr->tht >4.5)
		      {
			*return_code = CONIFERS_SUCCESS;
			
			      swo_hybrid_calc_d6_from_ht_and_dbh(return_code,
								 plant_ptr->tht,
								 plant_ptr->dbh,
								 &plant_ptr->d6,
								 c_ptr->d6_ht_dbh);
			    
			    if( *return_code != CONIFERS_SUCCESS)
			      {
				    plant_ptr->errors |= E_INVALID_D6;
				    error_count += 1;
			      }
		      }
		    else
		      {
		          /* if the dbh != 0 and tht < 4.5? */
			    *return_code = CONIFERS_SUCCESS;

			      swo_hybrid_calc_d6_from_total_height(  return_code, 
								     plant_ptr->tht, 
								     &plant_ptr->d6, 
								     c_ptr->d6_ht);
			    
			    if( *return_code != CONIFERS_SUCCESS)
			      {
				plant_ptr->errors |= E_INVALID_D6;
				error_count += 1;
			      }
		      }
		  }


          /* this is the code for computing the missing dbh when the tree */
          /* is taller than 4.5 feet                                        */
          /* if the plant record has a missing d12 and the total height is > 30 CM */
          /* this only applies to the CONIFERS_CIPS variant */
		if( plant_ptr->d12 == 0.0 && plant_ptr->tht > ( 30.0 * CM2FT ) )
		{
		    *return_code = CONIFERS_SUCCESS;
            
            //switch (variant)
			//{
			//    case CONIFERS_CIPS:
            //        plot_ptr = get_plot( plant_ptr->plot, n_points, plots_ptr );
            //        cips_calc_d12_from_ht_and_veg_cov(  return_code,
			//                    plant_ptr->tht,
			//                    plot_ptr->shrub_pct_cover,
			//                    &plant_ptr->d12,
			//                    c_ptr->dbh_ht_veg_cov );
            //        break;
			//	    
            //    default:
            //        break;
            //} 
        }

          /* this is the code for computing the missing dbh when the tree */
          /* is taller than 4.5 feet                                        */
		if( plant_ptr->d6 > 0.0 && plant_ptr->dbh == 0.0 && plant_ptr->tht > 4.5)
		  {
		    *return_code = CONIFERS_SUCCESS;
		    
			  swo_hybrid_calc_dbh_from_height_and_d6(return_code,
								 plant_ptr->d6,
								 plant_ptr->tht,
								 &plant_ptr->dbh,
								 c_ptr->d6_ht_dbh);
			
			if( *return_code != CONIFERS_SUCCESS)
			  {
			    plant_ptr->errors |= E_INVALID_DBH;
			    error_count += 1;
			  }
        }

		if(plant_ptr->expf <=0.0 && fixed_plot_radius > 0.0  )
		  {
		    if( plant_ptr->errors & E_INVALID_DBH)
		      {
			plant_ptr->expf = 0;
			plant_ptr->errors |= E_INVALID_EXPF;
			error_count += 1;
		      }
		    else if( plant_ptr->dbh > min_dbh )
		      {
			plant_ptr->expf = baf / (plant_ptr->dbh * plant_ptr->dbh * FC_I);
		      }
		    else
		      {
			plant_ptr->expf = SQ_FT_PER_ACRE / 
			  ( fixed_plot_radius * fixed_plot_radius * MY_PI );
		      }            

		    /* multiply the expansion factor by the number of   */
		    /* stems that this plant record represents          */
		    plant_ptr->expf *= plant_ptr->n_stems;
		  }
	      }

            /* if the plant is a shrub, then compute the d6 if it's missing */
	      if( is_shrub( c_ptr ) )
		{
		  if( plant_ptr->d6 == 0.0 )
		    {
		      /* compute the d6 from the height */
		      *return_code = CONIFERS_SUCCESS;
		      
			    swo_hybrid_calc_d6_from_total_height(  return_code, 
								   plant_ptr->tht, 
								   &plant_ptr->d6, 
								   c_ptr->d6_ht);
			  
			  if( *return_code != CONIFERS_SUCCESS)
			    {
			      plant_ptr->errors |= E_INVALID_D6;
			      error_count += 1; 
			    }
		    }

		}

	  }

	/* fill in the remaining values */
	plant_ptr->d6_area      = plant_ptr->d6 * plant_ptr->d6 * FC_I;
	plant_ptr->basal_area   = plant_ptr->dbh * plant_ptr->dbh * FC_I;
    plant_ptr->d12_area      = plant_ptr->d12 * plant_ptr->d12 * FC_I;
	    

	    /* now calculate the crown width for the plant record */
	    if( plant_ptr->crown_width <= 0.0 )
	      {
		*return_code = CONIFERS_SUCCESS;
		
		      swo_hybrid_calc_crown_width(   return_code,
						     plant_ptr->d6_area,
						     plant_ptr->tht,
						     &plant_ptr->crown_width,
						     &plant_ptr->crown_area,
						     c_ptr->crown_width );
		    
		    if( *return_code != CONIFERS_SUCCESS)
		      {
			plant_ptr->errors |= E_INVALID_CW;
			error_count += 1;
		      }
	      }

	    /*  MOD033  */
	    /* this should happen no matter what... */
	    plant_ptr->crown_area = plant_ptr->crown_width * plant_ptr->crown_width * MY_PI / 4.0;

	    /* fill in the expansion factors for those plant records    */
	    /* that don't have one filled in by using the sample        */
	    /* design records. Mostly, shrubs will have the expf filled */
	    /* by calling this function, but it should work for trees   */
	    if( plant_ptr->expf <= 0.0 )
	      {
		
		    swo_hybrid_calc_exp_from_cover_and_ca( return_code,
							   plant_ptr->pct_cover,
							   plant_ptr->crown_area,
							   &plant_ptr->expf );

            if( *return_code != CONIFERS_SUCCESS)
		    {
		      plant_ptr->errors |= E_INVALID_EXPF;
		      error_count += 1;
		    }

	      }
    
    }

  calc_plot_stats_2(    return_code, 
			n_species,
			species_ptr,
			n_coeffs,
			coeffs_ptr,
			n_plants,
			plants_ptr,
			n_points,
			plots_ptr );

    /* todo: need_error_trap_here */

    /* SECOND PASS */
    /* The second pass is required becuase the crown values are */
    /* dependent on the basic plot summary statistics which are */
    /* calculated for the plot before hand                      */
    /* now make the second pass */
    //error_count  = 0;
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
      {

	c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];


	/* don't go any further if you don't have a functional species code */
	if( c_ptr == NULL )
	  {
	    continue;
	  }

	/* set the pointer for the current plot that the plant is on */
	plot_ptr = get_plot( plant_ptr->plot, n_points, plots_ptr );
        
	/* fill in the percent cover for the plant recor */
	if( plant_ptr->pct_cover == 0.0)  /* fill in percent cover */
	  {
	    plant_ptr->pct_cover = 100.0 * plant_ptr->expf * plant_ptr->crown_area / SQ_FT_PER_ACRE;
	  }

	/* if the crown ratio is invalid or missing     */
	/* fill in the crown ratio for the plant record */
	if( is_tree( c_ptr ) )
	  {
	    if( plant_ptr->cr <= 0.0 || plant_ptr->cr > 1.0 )
	      {

		        //get_in_taller_attribs( plant_ptr, plot_ptr, bait, cait );
		        //cait_c       =   cait[CONIFER];
		        //cait_h       =   cait[HARDWOOD];
		        //cait_s       =   cait[SHRUB];
		        //*return_code = CONIFERS_SUCCESS;

		      swo_hybrid_calc_crown_ratio(return_code,  
						  plant_ptr->tht,
						  plant_ptr->d6,
						  &plant_ptr->cr,
						  c_ptr->crown_ratio);
	      }                              
	  }

	if( *return_code != CONIFERS_SUCCESS)
	  {
	    plant_ptr->errors |= E_INVALID_CR;
	    error_count += 1;
	  }

	plant_ptr->max_crown_width = 0.0;

	/*  MOD014 */
	  /* if it's a tree, compute the crown width */
	  if( is_tree( c_ptr ) )
	    {
	      *return_code = CONIFERS_SUCCESS;
	      
		    swo_hybrid_calc_max_crown_width(  return_code,
						      plant_ptr->dbh,
						      plant_ptr->tht,
						      &plant_ptr->max_crown_width,
						      c_ptr->max_crown_width);
     
		  if( *return_code != CONIFERS_SUCCESS)
		    {
		      plant_ptr->errors |= E_INVALID_MCW;   
		      error_count += 1;
		    }
	    }
      }
   
    if (error_count > 0)
      {
	*return_code = FILL_VALUES_ERROR;
	return;
      }
    *return_code = CONIFERS_SUCCESS;

}


/********************************************************************************/
/* swo_hybrid_project_plant                                                     */
/********************************************************************************/
/*  Description :   This function compute the plant growth for one year.        */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   September 12, 1999                                          */
/*  Returns     :   void                                                        */
/*  Comments    :   updated values are height_growth, d6_growth, dbh_growth,    */
/*                  crown_ratio_growth, crown_width_growth, and max_crown_width */
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
void swo_hybrid_project_plant(  
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   struct PLANT_RECORD     *plant_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           endemic_mortality,  
   int                     hcb_growth_on,      
   unsigned long           use_precip_in_hg,   
   unsigned long           use_rand_err )
{

   struct COEFFS_RECORD    *c_ptr;
   double                  awi;

   /* tree level stats */
   double  bat_total;
   double  bat_c;
   double  bat_h;
   double  bat_s;
   double  bat_c_h;
   double  cat_c;
   double  cat_h;
   double  cat_s;
   double  new_d6_area;
   double  new_d12_area;

/*    unsigned long htidx = 0; */

   double  normal;
   double  browse_random_unif_0_1;
   double  top_dam_random_unif_0_1;

   double  bait[PLANT_TYPES];
   double  cait[PLANT_TYPES];
    
   /* get the supporting structures for the plant */
   c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

   /* don't go any further if you don't    */
   /* have a functional species code       */
   if( c_ptr == NULL )
   {
      *return_code = CONIFERS_ERROR;
      return;
   }

//   Rprintf( "%s, %d\n", __FILE__, __LINE__ );

   /* get a uniform deviate for the browse */
   /* and one for the top damage           */
   normal                  = (double)gauss_dev();  
   browse_random_unif_0_1  = uniform_0_1();
   top_dam_random_unif_0_1 = uniform_0_1();

   /* calc the height growth...        */
   /* calc diam growth...              */
   /* calc crown recession...          */
   /* calc mortality...                */
   awi         =   plot_ptr->water_capacity * log( plot_ptr->mean_annual_precip );

    get_in_taller_attribs( plant_ptr, 
                            plot_ptr, 
                            bait, 
                            cait );

   bat_c       =   bait[CONIFER];
   bat_h       =   bait[HARDWOOD];
   bat_s       =   bait[SHRUB];

   cat_c       =   cait[CONIFER];
   cat_h       =   cait[HARDWOOD];
   cat_s       =   cait[SHRUB];

   bat_c_h     =   bat_c + bat_h;
   bat_total   =   bat_c + bat_h + bat_s;

   plant_ptr->tht_growth = 0.0;
   plant_ptr->d6_growth  = 0.0;
   plant_ptr->d12_growth  = 0.0;


   plant_ptr->cw_growth  = 0.0;
   plant_ptr->cr_growth  = 0.0;
   plant_ptr->dbh_growth = 0.0;
   plant_ptr->expf_change= 0.0;

	if( is_tree( c_ptr ) || is_shrub( c_ptr) )
	{

		/* coded jdh, july 31, 2011 */
		swo_hybrid_calc_height_growth( return_code,
			                plant_ptr->tht, 
			                plant_ptr->cr,
			                plant_ptr->d6_area,
			                
							/* the precip in the swohybrid model is *NOT* the	*/
							/* mean_annual_precip, but the growing season		*/
							/* precip, which includes only March, April,		*/
							/* and May.											*/
							//plot_ptr->mean_annual_precip, 
			                plot_ptr->growing_season_precip,
							
							plot_ptr->water_capacity,
			                plot_ptr->slope,
			                plot_ptr->aspect,
						
			                plot_ptr->ca_c,
			                cat_c,

			                plot_ptr->ca_h,         
			                plot_ptr->ca_s,         
            
			                cat_h,
			                cat_s,

			                plant_ptr->d6,
			                plot_ptr->elevation,
			                normal,                   
			                browse_random_unif_0_1,   
			                top_dam_random_unif_0_1,  
			                use_rand_err,         

			                species_ptr[plant_ptr->sp_idx].browse_damage,
			                species_ptr[plant_ptr->sp_idx].mechanical_damage,

			                use_precip_in_hg,     

							/* added arguments for the hybrid model */
			                species_ptr[plant_ptr->sp_idx].min_temp,
			                species_ptr[plant_ptr->sp_idx].max_temp,
			                species_ptr[plant_ptr->sp_idx].opt_temp,
							
							/* added plot level args for the radiation */
							plot_ptr->mean_monthly_temp,
			                plot_ptr->solar_radiation,

			                &plant_ptr->tht_growth,
			                c_ptr->ht_growth,
			                c_ptr->type); 

        if( *return_code != CONIFERS_SUCCESS )
        {

	        return;
        }

	}


	if( is_tree( c_ptr ) || is_shrub( c_ptr))
	{

		/* does this function have meaning here?				*/
		/* does it use the same equation as the swo variant???	*/
		swo_hybrid_calc_d6_growth(	return_code,
						plant_ptr->tht_growth,
						plant_ptr->crown_width,
						plant_ptr->tht,
						plant_ptr->d6,
                        cat_s,
						cat_c,
						cat_h,
						plot_ptr->ca_s,
						plot_ptr->ca_c,
						plot_ptr->ca_h,
						plot_ptr->water_capacity,
						&plant_ptr->d6_growth,
						c_ptr->d6_growth,
						c_ptr->type);

		if( *return_code != CONIFERS_SUCCESS )
        {
	         return;
        }

	}
   /* predict the new dbh from     */
   /* the new d6 and total height  */
   if( is_tree( c_ptr ) )
   {
	   swo_hybrid_calc_dbh_growth(  return_code,
									plant_ptr->tht,
									plant_ptr->tht_growth,
                                    plant_ptr->cr,
                                    plant_ptr->dbh,
                                    plant_ptr->d6,
                                    plant_ptr->d6_growth,
									&plant_ptr->dbh_growth,
									c_ptr->dbh_growth );

      if( *return_code != CONIFERS_SUCCESS )
      {
		return;
      }
            
      /* calculate the new crown ratio and    */
      /* calculate the difference             */
	  /* is this function valid anymore?		*/
      swo_hybrid_calc_cr_growth(   return_code,
			            hcb_growth_on,                
			            plant_ptr->tht,
			            plant_ptr->tht_growth,
			            plant_ptr->cr,
			            plot_ptr->ca_c,
			            plot_ptr->ca_h,
			            plot_ptr->ca_s,
			            uniform_0_1(),
			            &plant_ptr->cr_growth,
			            c_ptr->cr_growth);

      if( *return_code != CONIFERS_SUCCESS )
      {
	    return;
      }

      swo_hybrid_calc_max_crown_width( return_code,
			                plant_ptr->dbh + plant_ptr->dbh_growth,
			                plant_ptr->tht + plant_ptr->tht_growth,
			                &plant_ptr->max_crown_width,
			                c_ptr->max_crown_width);

      if( *return_code != CONIFERS_SUCCESS )
      {
	    return;
      }



   }

   /* calculate the new crown area for     */
   /* the plant record...  need to calc    */
   /* a temp new d6 area, total height,    */
   /* and crown width                      */
   new_d6_area = plant_ptr->d6 * plant_ptr->d6 * FC_I;
   new_d12_area = plant_ptr->d12 * plant_ptr->d12 * FC_I;
   
   
   if( is_tree( c_ptr ) || is_shrub( c_ptr))
   {

//  crown width growth hat
//  cwghat = (c0 + c2*(cw^0.5))*(hghat1)^c1

		swo_hybrid_calc_cw_growth( return_code,
				        plant_ptr->tht,
				        plant_ptr->tht_growth,
				        plant_ptr->crown_width,    
	                    plot_ptr->ca_c,
				        plot_ptr->ca_h,
					    plot_ptr->ca_s,
					    cat_c,
				        uniform_0_1(),
				        &plant_ptr->cw_growth,
					    c_ptr->cw_growth,
                        c_ptr->type);

        if( *return_code != CONIFERS_SUCCESS )
        {
	       return;
        }
   }

   if( endemic_mortality == 1 )    
   {

      swo_hybrid_calc_endemic_mortality(	return_code,
										  plant_ptr->expf,
										  &plant_ptr->expf_change,
										  &species_ptr[plant_ptr->sp_idx].endemic_mortality);

      if( *return_code != CONIFERS_SUCCESS )
      {
		return;
      }
   }


   *return_code = CONIFERS_SUCCESS;

}



/********************************************************************************/
/*                  calc_crown_width     S1                                     */
/********************************************************************************/
/*  Description :   calc_crown_width                                            */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -  return code for calling function to check    */
/*  d6_area                     -  basal area in inches of largest stem at base */
/*  total_height                -  total height of the plant                    */
/*  *pred_crown_area            -  predicted crown area                         */
/*  *pred_crown_width           -  predicted crown width in square feet         */
/*  vector<double> *coeffs_ptr  -  pointer to a vector of doubles that          */
/*                                   contain the coefficients for the           */
/*                                   functional species code                    */
/*                                                                              */
/********************************************************************************/
/*  Formula : log (ca) = b0 + b1 log (d6ba) + b2* log(h)   S1                   */ 
/*           where: ca      = crown area in square feet                         */
/*                  h       = height, feet                                      */
/*                  d6ba    = basal area of largest stem at 6 inches, in^2      */
/*  Source  : Uzoh and Ritchie 1996 PSW RP 227                                  */
/*  Coeffs  : CW   may need to add b3 due to minor LOF probs with df            */
/********************************************************************************/
void swo_hybrid_calc_crown_width( 
    unsigned long   *return_code,
    double          d6_area,
    double          total_height,
    double          *pred_crown_width,
    double          *pred_crown_area,
    double          *coeffs_ptr)
{

    double  b0;
    double  b1;
    double  b2;

    /* check for valid height */
    if(total_height <= 0.0)
    {
        *return_code = INVALID_INPUT_VAL;
        *pred_crown_width = 0.0;
        *pred_crown_area  = 0.0;
        return;
    }
    /* check for valid coefficients */
    if( coeffs_ptr == NULL )
    {
        /* MOD002   */
        *return_code = INVALID_COEFF;
        *pred_crown_width   = 0.0;
        *pred_crown_area    = 0.0;
        return;
    }

    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];
    b2 = coeffs_ptr[2];

    if( total_height < 0.51 )
    {
        *pred_crown_width = 0.25;
        *pred_crown_area  = 0.04908739;
    }
    else
    {

        *pred_crown_area = exp( b0 + 
                                b1 * log( d6_area * 144.0 ) + 
                                b2 * log( total_height ) );

		/* limit crown width to 60 feet */
		if(*pred_crown_area > 2827.0)
		{
			*pred_crown_area = 2827.0; 
		}

		*pred_crown_width = sqrt( *pred_crown_area * ONE_OVER_PI * 4.0); 
    }

    *return_code = CONIFERS_SUCCESS;

    if( *pred_crown_width < 0.0 )
    {
        *pred_crown_width   = 0.0;
        *pred_crown_area    = 0.0;
        *return_code        = CONIFERS_ERROR;
    }

}



/********************************************************************************/
/*                  calc_max_crown_width  S2                                    */
/********************************************************************************/
/*  Description :   calc_max_crown_width                                        */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  dbh                         -   dbh of the tree                             */
/*  total_height                -   total height of the subject tree            */
/*  *pred_max_crown_width       -   predicted MCW in feet                       */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : mcw = b0 + b1 * dbh + b2 * dbh * dbh                              */ 
/*           where: mcw= maximum crown width of an open grown tree in feet      */
/*                  dbh = breast height diameter in inches                      */
/*  Source  : Paine and Hann                                                    */
/*  Coeffs  : MW                                                                */
/********************************************************************************/

void swo_hybrid_calc_max_crown_width( 
    unsigned long   *return_code,
    double          dbh,
    double          total_height,
    double          *pred_max_crown_width,
    double          *coeffs_ptr    )
{
    double  b0;
    double  b1;
    double  b2;

    if( coeffs_ptr == NULL ) /* MOD002 */
    {
        *return_code = INVALID_COEFF;
        *pred_max_crown_width = 0.0;
        return;
    }
    
    if( dbh < 0.0 ) /* MOD013 */
    {
        *return_code = INVALID_INPUT_VAL;
        *pred_max_crown_width = 0.0;
        return;
    }

    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];
    b2 = coeffs_ptr[2];
    
    /* MOD012 */
    if( total_height <= 4.5 )
    {
        *pred_max_crown_width = b0 * (total_height / 4.5);
    }
    else
    {
        *pred_max_crown_width = b0 + b1 * dbh + b2 * dbh * dbh;
    }

    *return_code = CONIFERS_SUCCESS;

    if( *pred_max_crown_width < 0.0 )  /*  MOD006  */
    {
        *pred_max_crown_width   = 0.0;
        *return_code            = CONIFERS_ERROR;
    }
}



/********************************************************************************/
/*                  calc_crown_ratio (not for shrubs..)   S3                    */
/********************************************************************************/
/*  Description :   calc_crown_width                                            */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999  modified May 2008                         */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  d6                          -   basal diameter                              */
/*  double *pred_cr             -   pointer to predicted crown ratio            */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : CL = b0*H**b1 - exp( b2 + b3 d6/height)                           */ 
/*  Source  : Ritchie May 2008                                                  */
/*  Coeffs  : CR                                                                */
/********************************************************************************/
void swo_hybrid_calc_crown_ratio(
    unsigned long *return_code,  
    double  total_height,
    double  d6,
    double  *pred_cr,
    double  *coeffs_ptr    )
{

    double b0, b1, b3;
    double ratio;
    double newCR;           
    double crown_length;

    *return_code = CONIFERS_SUCCESS;

    if( coeffs_ptr == NULL )
    {
        *return_code = INVALID_COEFF;
        *pred_cr = 0.0;
        return;
    }

    if( total_height <= 0.0 || d6 <= 0.0)
    {
        *return_code = INVALID_INPUT_VAL;
        *pred_cr = 0.0;
        return;
    }

    ratio=d6/total_height;

    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];
    b3 = coeffs_ptr[3];

        crown_length = b0*pow(total_height,b1) * exp(b3*ratio);
    
    if (total_height <= 1.0)
    {
        crown_length = b0*total_height;
    }

    newCR  = crown_length/total_height;
    
    if( newCR < 0.0 )
    {
        newCR = 0.0;
    }
    if( newCR > 1.0 )
    {
        newCR = 1.0;
    }

    *pred_cr = newCR;

}



/********************************************************************************/
/*                  calc_d6_from_total_height       S4                          */
/********************************************************************************/
/*  Description :   Calculate a missing basal diameter from                     */
/*                  the total height.                                           */
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  *pred_d6                    -   predicted basal diameter from the function  */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula :   ln ( D6 ) = b0 * b1 * ln (H - 0.5)                              */
/*       Where   H=total height in feet                                         */
/*               D6=basal diameter at 6 inches above ground, in inches          */
/*                                                                              */
/*  Source: Ritchie, Static Height Diameter Equations S3                        */
/*  Coefficients: MS in coefficients file                                       */
/*  coefficients last fit in October 1999                                       */
/********************************************************************************/
void swo_hybrid_calc_d6_from_total_height(       
    unsigned long   *return_code,
    double          total_height, 
    double          *pred_d6,
    double          *coeffs_ptr    )
{
    double  b0;
    double  b1;

    /* initialize variables */
    *return_code    = CONIFERS_ERROR;
    *pred_d6        = 0.0;
    b0              = 0.0;
    b1              = 0.0;

    if( coeffs_ptr == NULL )
    {
        *pred_d6        = 0.0;
        *return_code    = INVALID_COEFF;
        return;
    }
        
    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];
    /*  MOD029  */
    if(total_height > 0.5)
    {
        *pred_d6 = exp( b0 + b1 * log( total_height - 0.5 ) );
        *return_code = CONIFERS_SUCCESS;
    }
    else
    {
        *pred_d6=0.1;
    }

    if( *pred_d6 < 0 )
    {
        *pred_d6    = 0.0; 
        *return_code        = CONIFERS_ERROR;
    }


}


/********************************************************************************/
/*                  calc_d6_from_ht_and_dbh          S5                         */
/********************************************************************************/
/*  Description :   Calculate a missing basal diameter from                     */
/*                  the total height and DBH.                                   */
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  *pred_d6                    -   predicted basal diameter from the function  */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula:  D6/dbh = 1.0/(b0 -b1 * exp ( -b2 ln(H - 4.5))                     */ 
/*       Where: D6=basal diameter at 6 inches measured in inches                */
/*              dbh= breast height diameter in inches                           */
/*              H = total tree or plant height of the largest stem in feet      */
/*  Source  : Ritchie Static Height Diameter Equations,   S5                    */
/*  Coeffs  : DH                                                                */
/********************************************************************************/
void swo_hybrid_calc_d6_from_ht_and_dbh(       
    unsigned long   *return_code,
    double          total_height,
    double          dbh,
    double          *pred_d6,
    double          *coeffs_ptr    )
{

    double  b0;
    double  b1;
    double  b2;

    /* perform error check for the correct number of coeffs */
    if( coeffs_ptr == NULL )
    {
        *return_code = INVALID_COEFF;
        return;
    }

    if( total_height <= 4.5 )
    {
        *return_code = INVALID_INPUT_VAL;
        return;
    }

    b0  = coeffs_ptr[0];
    b1  = coeffs_ptr[1];
    b2  = coeffs_ptr[2];

    /*   MOD014  */
    *pred_d6 = dbh / ( b0 - b1 * exp ( -b2 * ( total_height - 4.5 ) ) );
    *return_code = CONIFERS_SUCCESS;

    /* MOD006 */
    if( *pred_d6 < 0.0 )
    {
        *pred_d6    = 0.0;
        *return_code        = CONIFERS_ERROR;
    }
}


/********************************************************************************/
/*                  calc_dbh_from_total_height                   S6             */
/********************************************************************************/
/*  Description :   Calculate a missing basal diameter from total height        */
/*                  and dbh                                                     */
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   I'm not sure why this is in here! Inverted Larsen&Hann      */
/*  Arguments   :                                                               */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  *pred_dbh                   -   predicted basal diameter from the function  */
/*  vector<double> *coeffs_ptr  -   pointer to a vector of doubles that         */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula :  D = [(ln (H - 4.5) - b0) / b1]^(1 / b2)                          */ 
/*  Source  :  Larsen and Hann (1987), RP-49                                    */
/*  Coeffs  :  OD (params 0, 1, 2)                                              */
/********************************************************************************/
void swo_hybrid_calc_dbh_from_height(       
    unsigned long   *return_code,
    double          total_height,
    double          *pred_dbh,
    double          *coeffs_ptr    )
{

    double  b0;
    double  b1;
    double  b2;

    /*     MOD015     */
    if( coeffs_ptr == NULL || coeffs_ptr[2] == 0.0)
    {
        *return_code = INVALID_COEFF;
        return;
    }

    /*   MOD013       */
    if( total_height <= 4.5 )
    {
        *return_code = INVALID_INPUT_VAL;
        return;
    }

    b0  = coeffs_ptr[0];
    b1  = coeffs_ptr[1];
    b2  = 1.0/coeffs_ptr[2];
    /* This is inverted Larsen and Hann equation */

    *pred_dbh = pow( ( ( log( total_height - 4.5 ) - b0 ) / b1 ), b2 );
    *return_code = CONIFERS_SUCCESS;

    if( *pred_dbh < 0.0 )
    {
        *pred_dbh       = 0.0;
        *return_code    = CONIFERS_ERROR;
    }
}

/********************************************************************************/
/*                  calc_dbh_from_total_height_and_d6           S7              */
/********************************************************************************/
/*  Description :   Calculate a dbh from total height and basal diameter        */
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  d6                          -   basal diameter of the largest stem          */
/*  *pred_dbh                   -   predicted breast height diameter from the   */
/*                                  function                                    */
/*  vector<double> *coeffs_ptr  -   pointer to a vector of doubles that         */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula:  dbh = d6 * (b0 - b1 * exp (-b2* ln( H - 4.5)))                    */ 
/*  Source  : Ritchie, see d6_from_height_and_dbh                               */
/*  Coeffs  : DH                                                                */
/********************************************************************************/
void swo_hybrid_calc_dbh_from_height_and_d6(       
    unsigned long   *return_code,
    double          d6,
    double          total_height,
    double          *pred_dbh,
    double          *coeffs_ptr )
{

    double  b0;
    double  b1;
    double  b2;

    if( coeffs_ptr == NULL )
    {
        *return_code = INVALID_COEFF;
        return;
    }

    if( total_height <= 4.5 )
    {
        *return_code = INVALID_INPUT_VAL;
        return;
    }


    b0  = coeffs_ptr[0];
    b1  = coeffs_ptr[1];
    b2  = coeffs_ptr[2];
    
    *pred_dbh = d6 * ( b0 - b1 * exp ( -b2 * ( total_height - 4.5 ) ) );

    if( *pred_dbh < 0.0 )
    {
        *pred_dbh = 0.1;
    }

    *return_code = CONIFERS_SUCCESS;
}


/********************************************************************************/
/*                  calc_exp_from_cover_and_ca       S8                         */
/********************************************************************************/
/*  Description :   Equation to calculate expansion factor for shrubs (& hw)    */
/*                  entered on a plot+pct cover basis                           */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  pct_cover                   -   Percent cover on the plot for this observ.  */
/*  crown_area                  -   Crown area                                  */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/*  *pred_expf                  -   predicted expansion factor                  */
/********************************************************************************/
/*  Formula:  exp = (43560 * pct_cover / 100) / ca                              */ 
/*  Source  : Me                                                                */
/*  Coeffs  : none                                                              */
/********************************************************************************/
void swo_hybrid_calc_exp_from_cover_and_ca(
    unsigned long   *return_code,
    double          pct_cover,
    double          crown_area,
    double          *pred_expf)
{
    double  cover_total=0.0;
    *return_code        = CONIFERS_ERROR;
    
    /* Cover must be between 0 and 100 per record   */
    if( pct_cover < 0.0 || pct_cover > 100.0 )
    {
        *return_code    = CONIFERS_ERROR;   
        *pred_expf      = 0.0;
        return;
    }

    /* crown area must be > 0 to avoid bad things   */
    if( crown_area <= 0.0 )
    {
        *return_code    = CONIFERS_ERROR;   
        *pred_expf      = 0.0;
        return;
    }

    /* This is the total cover in sq ft from %cvr */
    cover_total     =   SQ_FT_PER_ACRE * pct_cover * 0.01;
    *pred_expf      =   cover_total / crown_area;           
    *return_code    =   CONIFERS_SUCCESS;
}



/********************************************************************************/
/*                  calc_height_from_d6             S9                          */
/********************************************************************************/
/*  Description :   Equation to total stem height from basal diamter            */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  unsigned long *return_code  -  return code for calling function to check    */
/*  d6                          -  basal diameter of largest stem               */
/*  *pred_total_height          -  predicted height diameter from the function  */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula:  H - 4.5 = exp ( ( log(d6) - b0) / b1 ) + 0.5                      */ 
/*  Source  :  Ritchie height diameter equations, inverted                      */
/*  Coeffs  :  MS                                                               */
/********************************************************************************/
void swo_hybrid_calc_height_from_d6(       
    unsigned long   *return_code,
    double          d6,
    double          *pred_total_height,
    double          *coeffs_ptr    )
{
    double  b0;
    double  b1;

    /* initialize variables */
    *return_code        = CONIFERS_ERROR;
    *pred_total_height  = 0.0;
    b0                  = 0.0;
    b1                  = 0.0;

    /*  MOD015   */
    if( coeffs_ptr == NULL || b1 == 0.0)
    {
        *pred_total_height  = 0.0;
        *return_code        = INVALID_COEFF;
        return;
    }

    /*   MOD013       */
    if( d6 <= 0.0 )
    {
        *pred_total_height  = 0.0;
        *return_code    = INVALID_INPUT_VAL;
        return;
    }

    /* set the values of the coeffs */
    b0  = coeffs_ptr[0];
    b1  = coeffs_ptr[1];
    
    *pred_total_height = exp( ( log( d6 ) - b0) / b1 ) + 0.5;
    *return_code = CONIFERS_SUCCESS;

    if( *pred_total_height < 0.0 )
    {
        *pred_total_height  = 0.0;
        *return_code        = CONIFERS_ERROR;
    }

}


/********************************************************************************/
/*                  calc_height_from_dbh          S10                           */
/********************************************************************************/
/*  Description :   Equation to total stem height from dbh                      */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  dbh                         -   dbh of the subject tree                     */
/*  *pred_total_height          -   predicted basal diameter from the function  */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/*                                                                              */
/********************************************************************************/
/*  Formula:      H = 4.5 + exp (b0 + b1 * dbh^b2)                              */ 
/*  Source  : Hann and Larsen (1987)                                            */
/*  Coeffs  : OD                                                                */
/********************************************************************************/
void swo_hybrid_calc_height_from_dbh(       
    unsigned long   *return_code,
    double          dbh,
    double          *pred_total_height,
    double          *coeffs_ptr    )
{

    double  b0;
    double  b1;
    double  b2;

    /* initialize variables */
    *return_code        = CONIFERS_ERROR;
    *pred_total_height  = 0.0;
    b0                  = 0.0;
    b1                  = 0.0;
    b2                  = 0.0;

    if( coeffs_ptr == NULL )
    {
        *pred_total_height = 0.0;
        *return_code = INVALID_COEFF;
        return;
    }

    /*   MOD013      */
    if( dbh <= 0.0)
    {
        *pred_total_height = 0.0;
        *return_code = INVALID_INPUT_VAL;
        return;
    }

    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];
    b2 = coeffs_ptr[2];

    *pred_total_height = 4.5 + exp( b0 + b1 * pow( dbh, b2 ) );
    *return_code = CONIFERS_SUCCESS;

    /* MOD006   */
    if( *pred_total_height < 0.0 )
    {
        *pred_total_height  = 0.0;
        *return_code        = CONIFERS_ERROR;
    }

}


/********************************************************************************/
/*           calc_num of stems per plant from cw_and_total_height     S11       */
/********************************************************************************/
/*  Description :   calc_nstems_from_cw_and_height                              */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   This function will round the predicted real                 */
/*                  value to a integer value                                    */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  crown_width                 -   crown width of the subject tree             */
/*  *pred_n_stems               -   predicted num of stems  from the function   */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/*                                                                              */
/*----------------------------------------------------------------------------- */
/*  Formula :  ln(#stems) = round (exp (b0 + b1 * cw + b2 cw*c2 +b3*height))    */ 
/*  Source  : Ritchie 11/99: Preliminary                                        */
/*  Coeffs  : ST                                                                */
/********************************************************************************/
void swo_hybrid_calc_nstems_from_cw_and_height(
    unsigned long   *return_code,
    double          total_height,
    double          crown_width,
    long            *pred_n_stems,
    double          *coeffs_ptr      )
{

    double  b0;
    double  b1;
    double  b2;
    double  b3;
    double  temp_stems;
    long    int_stems;

    /* initialize variables */
    *return_code        = CONIFERS_ERROR;
    *pred_n_stems       = 1;
        
    temp_stems          = 1.0;
    int_stems           = 1;

/*    MOD027    */
/*    general coefficients     */    
    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];
    b2 = coeffs_ptr[2];
    b3 = coeffs_ptr[3];


    /* temp_stems is a temporary value for holding the          */
    /* predicted missing n_stems the function is from           */
    /* a static n_stems predictive function STATIC: n_stems     */
    if( crown_width == 0.0 )
    {
        *return_code    =   CONIFERS_ERROR;
        return;
    }

    temp_stems    = exp(    b0 + 
                            b1 * crown_width + 
                            b2 * crown_width * crown_width + 
                            b3 * total_height   );
    
    /* force number of stems to be between 1 and 30     */
    if( temp_stems > 30.0 )
    {
        temp_stems = 30.0;
    }
    
    if( temp_stems < 1.0 )
    {
        temp_stems = 1.0;
    }

    /* convert value to a long integer from a double */
    int_stems       = (long) (temp_stems + 0.5);   
    *pred_n_stems   = int_stems;
    *return_code    = CONIFERS_SUCCESS;
}

/********************************************************************************/
/*           calc_volume        S12                                              */
/********************************************************************************/
/*  Description :   calc_vol4                                                   */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   January 11, 2000                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   This function will do volumes                               */
/*                                                                              */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  dbh                         -   crown width of the subject tree             */
/*  *volume4                    -   predicted volume to a 4 "                   */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/*                                                                              */
/*----------------------------------------------------------------------------- */
/*  Formula : volume = b0 * dbh ^ b1 * h ^b2 * dbh^b3                           */ 
/*  Source  : Wensel, Kirklely                                                  */
/*  Coeffs  : V4                                                                */
/********************************************************************************/
void swo_hybrid_calc_volume(
    unsigned long   *return_code,
    double          total_height,
    double          dbh,
    double          *pred_volume,
    double          *coeffs_ptr )
{
    double  b0;
    double  b1;
    double  b2;
    double  b3;
    
    /* initialize variables */
    *return_code        = CONIFERS_SUCCESS;
    *pred_volume        = 0.0;
    
	 if( coeffs_ptr == NULL )
    {
        *pred_volume = 0.0;
		*return_code = CONIFERS_ERROR;
        return;
    }
    

    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];
    b2 = coeffs_ptr[2];
    b3 = coeffs_ptr[3];

    if( b0 <= 0.0 || b1 <= 0.0 || b2 <= 0.0 )
    {
        *return_code    =   INVALID_INPUT_VAL;
        return;
    }

	if( dbh <= 0.0 || total_height <= 4.5)
	{
		*pred_volume = 0.0;
		return;
	}
	*pred_volume = b0 * 
			   pow(dbh, b1 ) * 
			   pow(total_height, b2) *
			   pow( b3, dbh);
    *return_code    = CONIFERS_SUCCESS;
	return;
}

/********************************************************************************/
/*           calc_biomass        S13                                            */
/********************************************************************************/
/*  Description :   calc_biomass                                                */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   January 11, 2000                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   This function will do biomass                               */
/*                  Found to be severely bogus in April 2005, Corrected by mwr  */
/*                   watch for corrections to english units in here, tricksy    */
/*                  since we are using multiple sources for eqations,           */
/*                  we compute all and then add. Since only one set of eqs      */
/*                  should be effective for any one species the others should   */
/*                  zero. We check for this in code below.                      */
/*                                                                              */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject plant           */
/*  basaldiameter               -   basal diameter of largest stem              */
/*  crown_width                 -   crown width of subject plant                */
/*  dbh                         -   breast height diameter                      */
/*  *biomass                    -   predicted biomass in tons                   */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/*                                                                              */
/*----------------------------------------------------------------------------- */
/*  Formula : Varies by species                                                 */
/*            biomass GG = sum(exp(a+b*Ln(dbh))*bk1)                            */
/*            biomass RW = exp(a+b*ln(cw)) *k2                                  */ 
/*            biomass TH = (b1 * ca +b2 )*bk1                                   */
/*  Source  : Gholz et al. (GG)                                                 */
/*             Ross and Walstad (RW)                                            */
/*             Tappeiner, Harrington & Walstad (TH)                             */
/*  Coeffs  : BM                                                                */
/********************************************************************************/
void swo_hybrid_calc_biomass(
    unsigned long   *return_code,
    double          total_height,
    double          basal_diameter,
	double          crown_width,
	double          dbh,
    double     		*pred_biomass,
    double          *coeffs_ptr )
{
    double  rs_b1;   /* R&W shrub intercept */
    double  rs_b2;   /* R&W shrub slope     */
    double  rs_b3;   /* R&W shrub mse for bias correction */

    double  rp_b1;   /* R&W pine intercept */
    double  rp_b2;   /* R&W pine slope     */
    double  rp_b3;   /* R&W pine mse for bias correction */


	double  gf_b1;   /* G&G Foliage intercept */
	double  gf_b2;   /* G&G Foliage slope */
	double  gb_b1;   /* G&G Branch intercept */
	double  gb_b2;   /* G&G Branch slope */
	double  gs_b1;   /* G&G Stem intercept */
	double  gs_b2;   /* G&G Stem slope */

	double  h_b1;    /* H&T Slope     */
	double  h_b2;    /* H&T intercept */

	double  r_bio;   /* biomass from Ross and Walstad  - add log bias correction */
	double  g_biof;  /* fol biomass from Gholz & Grier - already bias corrected  */
	double  g_biob;  /* branch mass from Gholz & Grier - already bias corrected  */
	double  g_bios;  /* stem   mass from Gholz & Grier - already bias corrected  */
	double  h_bio;   /* biomass from Harrington & Tappenier                      */
	
    
    /* initialize variables */
    *return_code        = CONIFERS_SUCCESS;
    *pred_biomass       = 0.0;
	r_bio = 0.0;
	g_biof = 0.0;
	g_bios = 0.0;
	g_biob = 0.0;
	h_bio  = 0.0;
	

	 if( coeffs_ptr == NULL )
    {
        *pred_biomass = 0.0;
		*return_code = CONIFERS_ERROR;
        return;
    }
	/* set coefficients    */
    rs_b1 = coeffs_ptr[0];
    rs_b2 = coeffs_ptr[1];
    rs_b3 = coeffs_ptr[2];

	rp_b1 = coeffs_ptr[3];
	rp_b2 = coeffs_ptr[4];
	rp_b3 = coeffs_ptr[5];

	gf_b1 = coeffs_ptr[6];
	gf_b2 = coeffs_ptr[7];
	gb_b1 = coeffs_ptr[8];
	gb_b2 = coeffs_ptr[9];
	gs_b1 = coeffs_ptr[10];
	gs_b2 = coeffs_ptr[11];
	
	h_b1  = coeffs_ptr[12];
	h_b2  = coeffs_ptr[13];

    if( dbh <= 0.0 && crown_width <=0.0 && total_height <= 0.0)
    {
		*pred_biomass = 0.0;
        *return_code  = INVALID_INPUT_VAL;
        return;
    }
	
	/* Ross & Walstad shrubs */
	if( rs_b1 < 0.0 && rs_b2 > 0.0 && rs_b3 > 0.0 && crown_width > 0.0)
	{
        r_bio = exp(rs_b1  + rs_b2 * log(crown_width * FT2CM) + rs_b3 * 0.5 ) * GRM2LB; 

		if(r_bio < 0.0)  /* this should NEVER happen */
		{
			*pred_biomass = 0.0;
			*return_code = CONIFERS_ERROR;
			return;
		}
		*pred_biomass = r_bio * LB2TON;
		*return_code  = CONIFERS_SUCCESS;
		return;
	}

	/* Ross & Walstad pine */
	if (rp_b1 < 0.0 && rp_b2 > 0.0 && rp_b3 > 0.0 && total_height > 0.0)
	{
		r_bio = exp(rp_b1  + rp_b2 * log(total_height * FT2CM) + rp_b3 * 0.5 ) * GRM2LB; 

		if( r_bio < 0.0 )  /* this should NEVER happen */
		{
			*pred_biomass = 0.0;
			*return_code = CONIFERS_ERROR;
			return;
		}
		*pred_biomass = r_bio * LB2TON;
		*return_code  = CONIFERS_SUCCESS;
		return;
	}

	/* Harrington and Tappeiner eqs for tanoak and madrone */
	if ( h_b1 > 0.0 && h_b2 < 0.0 && crown_width > 0.0)
	{
		h_bio= (h_b1 * (MY_PI*((crown_width*FT2M)*(crown_width*FT2M))/4.0) + h_b2) * KG2LB;
		if ( h_bio < 0.0 ) /* this can happen for small plants, set to arbitrary small value */
		{
			h_bio = 0.01;
		}
		*pred_biomass = h_bio * LB2TON;
		*return_code = CONIFERS_SUCCESS;
		return;
	}

	/* Gholz et al. biomass equations */
	if( gf_b1 < 0.0 && gb_b1 < 0.0 && gs_b1 < 0.0 && gf_b2 > 0.0 && gs_b2 && gb_b2 > 0.0)
	{
		if ( dbh <= 0.0)
		{
			*pred_biomass = 0.01 * LB2TON;
			*return_code = CONIFERS_SUCCESS;
			return;

		}
		
		g_biof	= exp(gf_b1 + gf_b2 * log( dbh * IN2CM )) * KG2LB;   /* Gholz Foliage Mass*/
		g_biob	= exp(gb_b1 + gb_b2 * log( dbh * IN2CM )) * KG2LB;   /* Gholz Branch Mass */
		g_bios	= exp(gs_b1 + gs_b2 * log( dbh * IN2CM )) * KG2LB;   /* Gholz Stem Mass   */
		
		*pred_biomass = (g_biof + g_biob + g_bios) * LB2TON;
		*return_code = CONIFERS_SUCCESS;
		return;
	}

	/* shouldn't be here unless there are no valid coefficients */
	*pred_biomass = 0.0;
	*return_code    = CONIFERS_ERROR;
	return;
}



/********************************************************************************/
/* Growth Models																*/
/********************************************************************************/


/********************************************************************************/
/*                  calc_d6_growth          D1                                  */
/********************************************************************************/
/*  Description :   calc_d6_growth                                              */    
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   March 17, 2004                                              */
/*  Returns     :   void                                                        */
/*  Comments    :   copied from swo.model by jeff hamann                        */
/*                  diameter growth for tree species...                         */
/*                  dg_trees=(b0*hg**2)*exp(b2*cw +b3*d6                        */
/*                            +b4*catcon +  b5*cathw + b6*catsh)                */
/*                                                                              */
/*                  diameter growth for shrub species...                        */
/*                  dg_shrubs= c0 + c1*hg + c2*height + c3*1/bat + c4*awi.      */
/*                                                                              */
/*                  for both combined                                           */
/*                  dg=dg_trees + dg_shrubs                                     */
/*  NOTE:     one of these two variables should equal 0, (all parms=0)          */
/*  Arguments   :                                                               */
/*  unsigned long *return_code - return code for calling function to check      */
/*  double  height_growth    - predicted height increment (annual)              */
/*  double  crown_width      - crown width (feet)                               */
/*  double  total_height     - total tree height (feet)                         */
/*  double  current_d6       - initial d6 (inches)                              */
/*  double  cat_shrubs       - crown area in taller shrubs sqftpr acre          */
/*  double  cat_conifers     - crown area in taller con sqftpr acre             */
/*  double  cat_hardwoods    - crown area in taller hw  sqftpr acre             */
/*  double  ca_shrubs        - crown area in shrubs sqft per acre               */
/*  double  ca_conifers      - crown area in conifer sqft per acre              */
/*  double  ca_hardwoods     - crown area in hw      sqft per acre              */
/*  double  *pred_d6_growth  - predicted dbh growth (inches)                    */
/*  vector<double> *coeffs_ptr - pointer to a vector of doubles that            */
/*                              contain the coefficients for the                */
/*                              functional species code                         */
/********************************************************************************/
/*  Formula : See above                                                         */ 
/*  Source  : Ritchie 03/04                                                     */
/*  Coeffs  : d6_growth                                                         */
/********************************************************************************/
void swo_hybrid_calc_d6_growth(
	unsigned long   *return_code,
    double          height_growth,
    double          crown_width,
    double          total_height,
    double          d6,
    double          cat_shrubs,
    double          cat_conifers,
    double          cat_hardwoods, 
    double          ca_shrubs,
    double          ca_conifers,
    double          ca_hardwoods,
	double          h20_holding_capacity,
    double          *pred_d6_growth,
    double          *coeffs_ptr,
    unsigned long   plant_type)
{

    double  dg_trees;
    double  dg_shrubs;

    double  b0;
    double  b1;
    double  b2;
    double  b3;
    double  b4;
    double  b5;
    double  b6;
    double  b7;
    double  b8;
    double  b9;
    double  b10;
    double  b11;
    double  b12;
    double  b13;
    double  b14;
    double  b15;

    double  c0;
    double  c1;
    double  c2;
    double  c3;
    double  c4;

    double  temp_cat_total;
    double  sum_cat;
    double  cat_hardwoods_conifers;
	double  whc;

    /* initialize variables */
    *return_code=CONIFERS_SUCCESS;
	dg_trees     = 0.0;
    dg_shrubs    = 0.0;
	whc = h20_holding_capacity;

    b0  = 0.0;
    b1  = 0.0;
    b2  = 0.0;
    b3  = 0.0;
    b4  = 0.0;
    b5  = 0.0;
    b6  = 0.0;
    b7  = 0.0;
    b8  = 0.0;
    b9  = 0.0;
    b10 = 0.0;
    b11 = 0.0;
    b12 = 0.0;
    b13 = 0.0;
    b14 = 0.0;
    b15 = 0.0;
    c0  = 0.0;
    c1  = 0.0;
    c2  = 0.0;
    c3  = 0.0;
    c4  = 0.0;


    sum_cat                = cat_shrubs + cat_conifers + cat_hardwoods; 
    cat_hardwoods_conifers = cat_conifers + cat_hardwoods;

    if( coeffs_ptr == NULL )
    {
        *pred_d6_growth = 0.0;
		*return_code=CONIFERS_ERROR;
        return;
    }

	if( crown_width <= 0.0 )
	{
		*pred_d6_growth =0.0;
		*return_code=INVALID_INPUT_VAL;
		return;
	}
/******************************if predicted height growht is negative set diameter growth to zero************************/
	if (height_growth <= 0.0)
	{
		*pred_d6_growth = 0.0;
		return;
	}


    if( sum_cat <= 0.0 )
    {
        temp_cat_total = 0.01;
    }
    else
    {
        temp_cat_total = sum_cat/SQ_FT_PER_ACRE;
    }

/**************************** this is now annual diameter growth for trees ***************************************/
	if( plant_type == CONIFER || plant_type == HARDWOOD)
	{
		
		/* these are from the original SWO model */
	    b0  = coeffs_ptr[0];
		b1  = coeffs_ptr[1];
		b2  = coeffs_ptr[2];
		b3  = coeffs_ptr[3];
		b4  = coeffs_ptr[4];
		b5  = coeffs_ptr[5];
		b6  = coeffs_ptr[6];
		b7  = coeffs_ptr[7];
		b8  = coeffs_ptr[8];
		b9  = coeffs_ptr[9];
		b10 = coeffs_ptr[10];
		b11 = coeffs_ptr[11];
		b12 = coeffs_ptr[12];
		b13 = coeffs_ptr[13];
		b14 = coeffs_ptr[14];
		b15 = coeffs_ptr[15];
/*  this is changed by mwr on  Sept 26, 2011 */
		dg_trees =   pow(height_growth, b1)
					*exp(   b0                           
							+	b2 * pow(crown_width, b6)      
							+	b3 * d6                     
							+	b4 * (cat_conifers/SQ_FT_PER_ACRE)*(cat_conifers/SQ_FT_PER_ACRE)    
							+	b5 * ( (cat_hardwoods/SQ_FT_PER_ACRE)*(cat_hardwoods/SQ_FT_PER_ACRE)  
							+	     + (cat_shrubs/SQ_FT_PER_ACRE)*(cat_shrubs/SQ_FT_PER_ACRE) )        
								);

		*pred_d6_growth = dg_trees;

	}
	else if(plant_type == SHRUB || plant_type == FORB)
	{
	    c0  = coeffs_ptr[0];
		c1  = coeffs_ptr[1];
		c2  = coeffs_ptr[2];
		c3  = coeffs_ptr[3];
		c4  = coeffs_ptr[4];

		dg_shrubs =    (  c0                              
						+ c1 * height_growth *2           
						+ c2 * total_height               
						+ c3 / (sqrt( temp_cat_total) )   
						+ c4 * whc) / 2.0;
		
		if( dg_shrubs < 0.0 )
		{
			dg_shrubs = 0.0;
		}
		
		*pred_d6_growth = dg_shrubs ;
	
	}
	else
	{
		*pred_d6_growth = dg_trees ;
	}

    if( *pred_d6_growth < 0.0 )
    {
        *pred_d6_growth = 0.0;
    }
}



/********************************************************************************/
/*                  swo_hybrid_calc_dbh_growth									*/
/********************************************************************************/
/*  Description :   swo_hybrid_calc_dbh_growth                                  */   
/*  Author      :   Martin W. Ritchie and Jeff D. Hamann						*/
/*  Date        :   July 31, 2011												*/
/*  Returns     :   void                                                        */
/*  Comments    :   reverted back to original Sept 2011                         */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double         total_height     - total tree height                         */
/*  double         height_growth    - predicted height increment (annual)       */
/*  double         crown_width		- crown width (feet)						*/
/*  double         catcon			- crown area in taller in conifers (ft^2)	*/
/*  double         cathwsh			- crown area in taller hardwds+shrubs (ft^2)*/

/*  double         min_temp			- min temp for phtotosynthasis				*/
/*  double         max_temp			- max temp for phtotosynthasis				*/
/*  double         opt_temp			- opt temp for phtotosynthasis				*/
/*  double         tday_c			- vector total number of growing days (above C)	*/
/*  double         srad				- vector solar radation						*/
/*  double         *pred_dbh_growth - predicted dbh growth (inches)             */
/*  vector<double> *coeffs_ptr      -   pointer to a vector of doubles that     */
/*                                      contain the coefficients for the        */
/*                                      functional species code                 */
/********************************************************************************/
/*  Formula :																	*/ 
/*	     dghat  = exp(b0														*/
/*             + b1 * log(hghat1)												*/
/*             + b2 * cw^0.5													*/
/*             + b3 * basal_diameter (d6)										*/
/*             + b4 * (catcon)^2												*/
/*             + b5 * (cathwsh)^2												*/
/*             + b6 * log(SRT))													*/
/*  Source  :                                                                   */
/*  Coeffs  :																	*/
/********************************************************************************/
void swo_hybrid_calc_dbh_growth(   
    unsigned long   *return_code,
    double          total_height,
    double          height_growth,
    double          crown_ratio,
    double          current_dbh,
    double          current_d6,
    double          d6_growth,
    double          *pred_dbh_growth,
    double          *coeffs_ptr )
{

    double  temp_dbh;
    double  b0;
    double  b1;
    double  b2;

	*return_code=CONIFERS_SUCCESS;

    if( coeffs_ptr == NULL )
    {
        *pred_dbh_growth    = 0.0;
        *return_code        = INVALID_COEFF;
        return;
    }

    if( d6_growth < 0.0)
    {
        *pred_dbh_growth    = 0.0;
        *return_code        = INVALID_INPUT_VAL;
        return;
    }

    temp_dbh    = 0.0;
    
    b0          = coeffs_ptr[0];
    b1          = coeffs_ptr[1];
    b2          = coeffs_ptr[2];
    
    if( total_height > 4.5 ) 
    {
        temp_dbh  = d6_growth * ( b0 + exp( b1 + b2 * total_height));
    }
    else
    {
        temp_dbh = 0.0;
    }

    *pred_dbh_growth    = temp_dbh;
    *return_code        = CONIFERS_SUCCESS;

    if( *pred_dbh_growth < 0 )
    {
        *pred_dbh_growth    = 0.0;
        *return_code        = CONIFERS_ERROR;
    }
}



/********************************************************************************/
/*                  calc_height_growth        D3                                */
/********************************************************************************/
/*  Description :   calc_height_growth                                          */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   March 17, 2004                                              */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  crown_ratio                 -   crown ratio of the subject tree             */
/*  d6_area                     -   basal diam (sqft) (at 6") of lrgest stem    */
/*  precip                      -   annual precip in inches                     */
/*  h20_holding_capacity        -   water holding capacity (inches)             */
/*  slope                       -   slope of the site in percent                */
/*  aspect                      -   aspect of the site in degrees               */
/*  hardwood_ba                 -   hardwood basal area largest stem @ 6 inches */
/*  hardwood_ca                 -   hardwood crown area sqft per acre           */
/*  shrub_ba                    -   shrub basal area largest stem @ 6 inches    */
/*  shrub_ca                    -   shrub crown area sq ft per acre             */
/*  conifer_ba                  -   conifer basal area at 6 inches sqft pr acre */
/*  conifer_ca                  -   conifer crown area sqft per acre            */
/*  bainl_c                     -   basal area in taller conifers at 6 inches   */
/*  bainl_hw                    -   basal area in taller hw at 6 inches         */
/*  bainl_sh                    -   basal area in taller shrubs at 6 inches     */
/*  basal_d,                    -   basal diameter at 6", in inches             */
/*  elevation                   -   elevation in feet                           */
/*  random_norm_0_1             -   normal random variable with mean 0, var 1   */ 
/*  random_unif_0_1             -   uniform random var. with range 0, 1         */
/*  random_unif_0_1a            -   a second unif ran. var with range 0, 1      */
/*  ind_random                  -   indicator for random error (0 or 1) 1=yes   */
/*  prob_browse                 -   probability of browse, only if h< 4.5       */
/*  prob_top_damage             -   probability of top damage                   */
/* -----------------------------------------------------------------------------*/
/*  *height_growth              -   predicted height growth                     */
/*  vector<double> *coeffs_ptr  -   pointer to a vector of doubles that         */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : ln(hg)= XB   for trees                                            */
/*               hg = XB   for shrubs                                           */
/*  Source  : Ritchie Fit 11/1999; 04/2000; 12/2001; 08/2006                    */
/*  Coeffs  : HG                                                                */
/********************************************************************************/
void swo_hybrid_calc_height_growth(
	unsigned long   *return_code,	
    double          total_height, 
    double          crown_ratio,
    double          d6_area,
    double          precip, 
    double          h20_holding_capacity,
    double          slope,
    double          aspect,
    double          cacon,
    double          catcon,
    double          cahw,
    double          cash,
    double          cathw,
	double          catsh,
	double          basal_d,
    double          elevation,
    double          random_norm_0_1,
    double          random_unif_0_1a,
    double          random_unif_0_1b,
    long            ind_random,
    double          prob_browse,
    double          prob_top_damage,
    unsigned long   use_precip, 
    
	double			min_temp,		/* species specific minumum temperature, in C	*/
	double			max_temp,		/* species specific maximum temperature, in C	*/
	double			opt_temp,		/* species specific optimal temperature, in C	*/
	double			*tday_c,		/* tday_c is a vector of monthly temps, in C	*/
	double			*srad,			/* solar radiation */

    double          *height_growth,
    double          *coeffs_ptr,
	unsigned long   plant_type)
{

/****************************  right now b25 and b26 are spares ******************************************/
    double  a0,  a1,  a2,  a3,  a4,  a5,  a6,  a7,  a8,  a9, a10;

	
	/* these are for the browse and damage? */
	double				b10, b11;

    double  whc;
	double  height_var;
	double  height_for_error;
    int     broken;
    int     browsing;
    double  cat;
    double  temp_growth;

	/* temp variables for development */
	//double	monthly_temp = 20.0;
	double	power_term;				/* another temp number	*/
	double	t1num;					/* temp number?			*/
	double	t2num;					/* temp number?			*/
	//double	tday_c;					/* which is what?		*/
	double	gspcp;					/* which is what?		*/
	double	SRT;					/* which is what?		*/
	double	tmod[12];					/* which is what?		*/
	
	double	hcb1;

	unsigned long	i;

/**************************** initialize variables *******************************************************/
	*return_code=CONIFERS_SUCCESS;

    if( coeffs_ptr == NULL )
    {
        *height_growth = 0.0;
		*return_code   = CONIFERS_ERROR;
        return;
    }

/*************************** check for some errors *******************************************************/
    if( total_height < 0.0)
    {
        *height_growth    = 0.0;
        *return_code      = INVALID_INPUT_VAL;
        return;
    }

    if(total_height <= 0.5)                                       
    {
        *height_growth = 0.20;
		return;
    }

/*********************************************************************************************************/
	*height_growth   = 0.0;                               /* set to zero for default                     */
    temp_growth      = 0.0;

	whc              = h20_holding_capacity;
	   	
	height_var       = random_norm_0_1;                   /* height growth error by default is N(0,1)    */

	height_for_error = total_height;

    a0                  = 0.0;
    a1                  = 0.0;
    a2                  = 0.0;
    a3                  = 0.0;
    a4                  = 0.0;
    a5                  = 0.0;
    a6                  = 0.0;
    a7                  = 0.0;
    a8                  = 0.0;
    a9                  = 0.0;
    a10                 = 0.0;
    

    broken              = 0;
    browsing            = 0;

	SRT = 0.0;

	memset( tmod, 0, sizeof( double ) * 12 );
	

    cat = (catcon+cathw+catsh)/SQ_FT_PER_ACRE;


/****************************** set up the browse coefficient (indicator) ***********************************/
    if(prob_browse > random_unif_0_1a && total_height <= 4.5)
    {
        browsing = 1;
    }

/******************************  set up the top damage coefficient (indicator) ******************************/
    if(prob_top_damage > random_unif_0_1b)
    {
        broken   = 1;
    }

/******************************  turn off the random error component if ind_random=0 ************************/
    if(ind_random == 0) 
	{
		height_var = 0.0;   
    }

/****************************** if broken then shut off random error ****************************************/
	if (broken) 
	{
		height_var = 0.0;
	}
	
/****************************** if browsed shut off random error ********************************************/
	if (browsing) 
	{
		height_var = 0.0;
	}

/*********************  the next stuff will keep random error acting ok. ***********************************/
/*********************	first we will make sure that we dont get a huge error term, then we  ***************/
/*********************	will limit the random normal component to between the 10th and 90th percentile******/
    if(total_height > 15.0 )
	{
		height_for_error = 15.0;
	}
	if(height_var < -1.645)                                         /*  keep random normal within reason     */
	{
		height_var = 0.0;                                        /*  if large negative, set to .0       */
	}
	if(height_var > 1.645)
	{
		height_var = 0.0;                                         /*   if large positive, set to .0      */
	}

/**********************  make sure that log of crown ratio doesn't go goofy **********************/
    if ( crown_ratio < 0.01 )
    {
        crown_ratio = 0.01;
    }

    if( plant_type == CONIFER || plant_type == HARDWOOD ) 
    {
				
		//## parameter  for DF
		//a0=  -1.79173;
		//a1=   0.041419;
		//a2=  0.736826;
		//a3=  -0.02037;
		//a4=   0.533997;
		//a5=  -0.17748;
		//a6=  -0.16519;
		//a8=   0.014363;
		//a9=   0.230168;

		a0          = coeffs_ptr[0];
		a1          = coeffs_ptr[1];
		a2          = coeffs_ptr[2];
		a3          = coeffs_ptr[3];
		a4          = coeffs_ptr[4];
		a5          = coeffs_ptr[5];
		a6          = coeffs_ptr[6];
		a7          = coeffs_ptr[7];
		a8          = coeffs_ptr[8];
		a9          = coeffs_ptr[9];
        a10         = coeffs_ptr[10];
	/* for now, just comment these out */
	b10 = (double)coeffs_ptr[11] * (double)broken;
	b11 = (double)coeffs_ptr[12] * (double)browsing;


	/* this is a species level variable */
	power_term = ( max_temp - opt_temp ) / ( opt_temp - min_temp );

	/* loop over the months to compute the solar radiation temperature modifier. */
	for( i = 0; i < 12; i++ )
	    {

		/* tday_c is a vector */
		t1num = tday_c[i] - min_temp;
		if( t1num < 0.0 )
		    {
			    t1num = 0.0;
		    }

		// t2numdf=tmaxdf - tday_c;  if t2numdf<0 then t2numdf=0;
		t2num = max_temp - tday_c[i];
		if( t2num < 0.0 )
		    {
			    t2num = 0.0;
		    }

		/* this is a vector??? */
		    tmod[i] = pow( ( t1num / ( opt_temp - min_temp ) ) * ( t2num / ( max_temp - opt_temp ) ), power_term );

		    SRT += srad[i] * tmod[i] / 1000.0;

	    }

	/* gspcp ?= growing_season_precip ?= precip				*/
	/* input is in inches/year, model version in mm/year	*/	
	gspcp = precip; 

	hcb1 = total_height * ( 1.0 - crown_ratio );	/* units */

	/* compute the height growth, in feet. */
	/* this test returns 437.827, which is probably incorrect.		*/
	/* this test returns 1.3260284159948714, which is more like it	*/
    /* changed the crown area to unitless 9/2011 & updated function mwr*/
	temp_growth = exp(a0
			        + a1 * whc
			        + a2 * log(total_height)
			        + a3 * ( pow( total_height, a10 ) ) * 0.1
			        + a4 * log( crown_ratio )
			        + a5 * (catcon/SQ_FT_PER_ACRE) * (catcon/SQ_FT_PER_ACRE)
			        + a6 * ( cahw/SQ_FT_PER_ACRE )*  (cahw/SQ_FT_PER_ACRE)
                    + a7 * ( cash/SQ_FT_PER_ACRE )*( cash/SQ_FT_PER_ACRE )
			        + a8 * log(gspcp)
			        + a9 * ( SRT ) );


        if( (temp_growth) <= 0.0 )
        {
            *height_growth  = 0.0 + (b10 + b11);
        }
        else
        {
	        *height_growth  = temp_growth + (b10 + b11);          /* damage adjustments                    */
        }
    }
    else if(plant_type == SHRUB )        
    {
		a0  = coeffs_ptr[0];
		a1  = coeffs_ptr[1];
		a2  = coeffs_ptr[2];
		a3  = coeffs_ptr[3];
		a4  = coeffs_ptr[4];

		*height_growth= a0/total_height 
                  + exp(a1 + a2*(log(total_height)) + a3*total_height*basal_d + a4*cat*cat);

        /****for unknown brush species set growth to zero *****/
        if (a1 == 0 && a2 == 0 && a3 == 0 && a4 == 0)
        {
            *height_growth=0.0;
        }
    }
	else
	{
		*height_growth=0.0;
	}

	/******************** this concludes the tree height growth model  ******************************************/
    
	/* if predicted height is <0.5 then        */
	if( total_height + *height_growth <= 0.5 )                       
    {
        *height_growth = 0.0 ;
    }
    
	return;

}



/********************************************************************************/
/*                  calc_cr_growth                                              */
/********************************************************************************/
/*  Description :   change in crown ratio                                       */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   April 21, 2000                                              */
/*  Returns     :   void                                                        */
/*  Comments    :   two stages (1) prob of recession for 2 years                */
/*                             (2) change in hcb for 2 years                    */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long  *return_code     - return code for calling function to check */
/*  int            cr_flag          - flag to turn on hcb 0=no 1=yes            */
/*  double         total_height     - total tree height                         */
/*  double         height_growth    - predicted height increment (annual)       */
/*  double         crown_ratio      - crown ratio                               */
/*  double         conifer_ca       - conifer crown area                        */
/*  double         hardwood_ca      - hardwood crown area                       */
/*  double         shrub_ca         - shrub crown area                          */
/*  double         unif_0_1         - uniform random variate between 0 and 1    */
/*  double         *cr_growth       - pointer to the change in cr               */
/*  vector<double> *coeffs_ptr      - pointer to a vector of doubles that       */
/*                                    contain the coefficients for the          */
/*                                    functional species code                   */
/********************************************************************************/
/*  Formula : hcbgro = b0 + b1 * cl + b2 * conca + b3 hwca + b4 * shca          */
/*            prob_hgro = exp( bx ) / ( 1 + exp(bx) )                           */ 
/*  Source  : Ritchie April 2000                                                */
/*  Coeffs  : CG                                                                */
/********************************************************************************/
void swo_hybrid_calc_cr_growth(
    unsigned long   *return_code,
    int             hcb_growth_on,
    double          total_height,
    double          height_growth,
    double          crown_ratio,
    double          conifer_ca,
    double          hardwood_ca,
    double          shrub_ca,
    double          uniform_0_1,
    double          *cr_growth,
    double          *coeffs_ptr )
{
    double      b0;              /* temporary coefficients                   */
    double      b1;
    double      b2;
    double      b3;
    double      b4; 
    double      b5;
    double      b6;
    double      b7;
    double      b8;
    double      b9;
    double      hcb_growth;      /* predicted hcb grow if crown receeds      */
    double      crown_length;    /* starting crown length                    */
    double      temp_exponent;   /* a temporary variable for logistic func   */ 
    double      prob_hcb;        /* probability hcb does not receed          */

/*   set the temporary variables     */

    hcb_growth  = 0.0;
    crown_length= crown_ratio * total_height;

/*    set the coefficients           */

    b0          = coeffs_ptr[0];
    b1          = coeffs_ptr[1];
    b2          = coeffs_ptr[2];
    b3          = coeffs_ptr[3];
    b4          = coeffs_ptr[4];
    b5          = coeffs_ptr[5];
    b6          = coeffs_ptr[6];
    b7          = coeffs_ptr[7];
    b8          = coeffs_ptr[8];
    b9          = coeffs_ptr[9];


    if( coeffs_ptr == NULL )
    {
        *cr_growth    = 0.0;
        *return_code  = INVALID_COEFF;
        return;
    }

    if( !hcb_growth_on )  /*  on - off switch because hcb change is biannual */
                          /*  off = no change in hcb, i.e. change in hcb=0   */
    {
        *cr_growth    = ( crown_length + height_growth ) / 
                        ( total_height + height_growth ) -
                        ( crown_length / total_height  );
        *return_code  = CONIFERS_SUCCESS;
        return;
    }


/* (1)   first predict change in height to crown base              */    

    hcb_growth  = b0 
                + b1 * crown_length 
                + b2 * 0.00001 * conifer_ca
                + b3 * 0.00001 * hardwood_ca
                + b4 * 0.00001 * shrub_ca;

    if( hcb_growth < 0.0 )          /*  this should NEVER happen!  */
    {
        hcb_growth   = 0.0;
        *cr_growth   = 0.0;
        *return_code = CONIFERS_ERROR;
		return;
    }
    if( hcb_growth > ( crown_length + height_growth ))
    {
    /*  MOD039   added the 0.1 fudge factor to keep cr non-zero & positive */
        hcb_growth   = crown_length + height_growth - 0.1; /*  may never happen! */
    }

/* (2)   next calculate the probability of change in crown base    */

    temp_exponent   = exp( b5 
                         + b6 * crown_ratio  
                         + b7 * conifer_ca  / 43560.0 
                         + b8 * hardwood_ca / 43560.0
                         + b9 * shrub_ca    / 43560.0);

/* (2a)  logistic function                                         */

    prob_hcb        = temp_exponent / ( 1.0 + temp_exponent );

/* (3)   if the uniform r.v. is greater then apply change in cr    */

    if( uniform_0_1 > prob_hcb )  
    {
        *cr_growth   = ( crown_length + height_growth   - hcb_growth ) / 
                       ( total_height + height_growth ) -
                       ( crown_length / total_height  );

        *return_code = CONIFERS_SUCCESS;
    }
    else
    {
        *cr_growth   = ( crown_length + height_growth ) / 
                       ( total_height + height_growth ) -
                       ( crown_length / total_height  );
/*      *cr_growth   = 0.0;   */
        *return_code = CONIFERS_SUCCESS;
    }
}


/********************************************************************************/
/*                  calc_cw_growth         D6                                   */
/********************************************************************************/
/*  Description :   calc_cw_growth                                              */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   March 17, 2004                                              */
/*  Returns     :   void                                                        */
/*  Comments    :   This is a new function                                      */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double         height_growth    - predicted height increment (annual)       */
/*  double         crown_width      - crown width (feet)                        */
/*  double         ca_conifers      - crown area in conifers  (ft2/acre)        */
/*  double         ca_hardwoods     - crown area in hardwoods (ft2/acre)        */
/*  double         ca_hardwoods     - crown area in shrubs    (ft2/acre)        */
/*  double         *pred_cw_growth  - predicted dbh growth (feet)               */
/*  vector<double> *coeffs_ptr      -   pointer to a vector of doubles that     */
/*                                      contain the coefficients for the        */
/*                                      functional species code                 */
/********************************************************************************/
/*  Formula : cwgro = (c0 + c2 * root cw + c3 cacon + c4 cahw2+c5*cash2)*hg**c1 */ 
/*  Formula : cwgro = hgrow * (c8*height^c9)* exp( c10*catcon)                  */
/*  Source  : Ritchie 11/99                                                     */
/*  Coeffs  : CG                                                                */
/********************************************************************************/
void swo_hybrid_calc_cw_growth(   
	unsigned long   *return_code,
	double          total_height,
	double		    height_growth,
	double		    crown_width,
	double		    ca_conifers,
	double		    ca_hardwoods,
	double		    ca_shrubs,
	double		    catcon,
	double		    unif_0_1,
	double          *pred_cw_growth,
	double          *coeffs_ptr,
    unsigned long   plant_type)
{

	double  temp_cwg;
	double  c0;
	double  c1;
	double  c2;
	//double  c3;
	//double  c4;
	//double  c5;
	//double  c6;
	//double  c7;
	//double  c8;
	//double  c9;
	//double  c10;

	*return_code	= CONIFERS_SUCCESS;

    if( coeffs_ptr == NULL )
    {
        *pred_cw_growth    = 0.0;
        *return_code        = INVALID_COEFF;
        return;
    }
	if(ca_conifers < 0.0 || ca_hardwoods < 0.0 || ca_shrubs < 0.0)
	{
		*pred_cw_growth	=0.0;
		*return_code	= INVALID_INPUT_VAL;
		return;
	}

	if(crown_width <= 0.0 )
	{
		*pred_cw_growth	=0.0;
		*return_code	= INVALID_INPUT_VAL;
		return;
	}

	if( height_growth < 0.0 ) 
    {
        *pred_cw_growth = 0.0;
		*return_code	= CONIFERS_SUCCESS;
		return;
    }


    temp_cwg    = 0.0;
    
	//c0=   0.836485;
	//c1=   0.944264;
	//c2=  -0.15986;

	c0      = coeffs_ptr[0];
	c1      = coeffs_ptr[1];
	c2      = coeffs_ptr[2];

	/*  for trees,  hw & conifers c8 is an indicator parameter 0 for tree */
	if( plant_type == CONIFER || plant_type == HARDWOOD )
	{

		temp_cwg = ( c0 + c2 * sqrt( crown_width ) ) * pow( height_growth, c1 );
 
		if( temp_cwg < 0.0 ) /* if pred growth is negative */
		{
			*pred_cw_growth    = 0.00;  /* Growth will be set to 0 */
			*return_code        = CONIFERS_SUCCESS;
			return;
		}

		*pred_cw_growth     = temp_cwg ;
		*return_code        = CONIFERS_SUCCESS;
		return;
	}


	//if(plant_type == SHRUB)   /* if it is a shrub, change to a different function using c8 c9 c10*/
	//{
	//    temp_cwg = (height_growth*2.0) * c8*pow(total_height, c9) * exp(c10*catcon);
    //    if(temp_cwg < 0.0)
    //    {
    //        *pred_cw_growth =0.0;
    //        *return_code    = CONIFERS_SUCCESS;
    //        return;
    //    }
	//    *pred_cw_growth  = 0.5* temp_cwg;
	//    *return_code     = CONIFERS_SUCCESS;
	//    return;
	//}
    //else
   // {
   //     *pred_cw_growth  = 0.0;
//	    *return_code     = CONIFERS_SUCCESS;
//	    return;
//    }


}



/********************************************************************************/
/*                  calc_endemic_mortality                                      */
/********************************************************************************/
/*  Description :   calculate the background mortality                          */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   January 3, 2000                                             */
/*  Returns     :   void                                                        */
/*  Comments    :   This should be called only if user indicates it is ON       */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double        expansion factor  - current expansion factor                  */
/*  double        *pred_mortality   - predicted mortality                       */
/*  vector<double> *coeffs_ptr      -   pointer to a vector of doubles that     */
/*                                      contain the coefficients for the        */
/*                                      functional species code                 */
/********************************************************************************/
/*  Formula : expansion = expansion*(1- mortality rate)                         */ 
/*  Source  : SYSTUM-1                                                          */
/*  Coeffs  : MO                                                                */
/********************************************************************************/
void swo_hybrid_calc_endemic_mortality(   
    unsigned long   *return_code,
    double          expansion_factor,
    double          *pred_mortality,
    double          *coeffs_ptr )
{

    double  b0;
    double  current_mort;
    /*  initialize parameters    */
    /* *pred_mortality = 0.0;    */
    b0              = coeffs_ptr[0];
    current_mort	= *pred_mortality;

    if( coeffs_ptr == NULL )
    {
        *pred_mortality     = 0.0;
        *return_code        = INVALID_COEFF;
        return;
    }

    /*   check that the starting expansion factor and mortality rate are in bounds  */
    if( expansion_factor < 0.0 || b0 < 0.0 || b0 > 0.10 || expansion_factor < current_mort)
    {
        *pred_mortality     = 0.0;
        *return_code        = INVALID_INPUT_VAL;
        return;
    }
    
    /*   this will calculate the number of trees to die using the endemic mortality and current mort if any */
    *pred_mortality     = (expansion_factor - current_mort)* b0;
    *return_code        = CONIFERS_SUCCESS;

    if( *pred_mortality < 0.0 )
    {
        *pred_mortality     = 0.0;
        *return_code        = CONIFERS_ERROR;
    }

}
