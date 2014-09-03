/********************************************************************************/
/*                                                                              */
/*  stats.c                                                                     */
/*  functions used to predict values for the CONIFERS growth model              */
/*                                                                              */
/********************************************************************************/
                                                                                
/*  $Id: stats.c 856 2012-02-21 23:53:56Z hamannj $  */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>          /*  this is only here for testing   */

#include "conifers.h"

static int compare_summaries_by_code( 
				     const void *ptr1, 
				     const void *ptr2 );

static int compare_trees_by_plot(
				 const void *ptr1,
				 const void *ptr2 );


static int compare_htn_by_plant_tht_expf( 
					 const void *ptr1, 
					 const void *ptr2 );

struct HTN_RECORD {
  int is_tree;
  double    tht;
  double    expf;
};

/* todo: these need to go into the generic file/function for all variants */
static   void calc_volume(
      unsigned long   *return_code,
      double          total_height,
      double          dbh,
      double          *pred_volume,
      double          *coeffs_ptr);

static   void calc_biomass(
      unsigned long   *return_code,
      double          total_height,
      double          basal_diameter,
      double          crown_width,
      double		  dbh,
      double     	  *pred_biomass,
      double          *coeffs_ptr );



/* put the get.set attribs_in_taller functions here */

/********************************************************************************/
/* fill_in_whc_and_precip                                                       */
/********************************************************************************/
/*  Description :   This function fills in the missing values for the           */
/*                  sampled stand. The primary function is to convert from      */
/*                  SI to precip and whc.                                       */
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   February  04, 2000                                          */
/*  Returns     :   void                                                        */
/*  Comments    :   Modified in December 2001 to fix problems with site         */
/*                  conversions with whc                                        */
/*                  ANY CHANGES TO THIS FUNCTION MAY REQUIRE CHANGES TO THE     */
/*                  calc_sites_from_whc function to match                       */
/*  Arguments   :   unsigned long *return_code  - pointer to a return code      */
/*                  struct SAMPLE_RECORD *sample_ptr                            */
/*                  double dfsi    - douglas-fir site index                     */
/*                  double ppsi    - ponderosa pine                             */
/*                  double *whc    - pointer to water holding capacity          */
/*                  double *precip - pointer to annual precipitation            */
/********************************************************************************/
void fill_in_whc_and_precip( 
			    unsigned long           *return_code,
			    double                  dfsi,
			    double                  ppsi,
			    double                  *whc,
			    double                  *precip)
{

  double temp_whc;

   
  /* force these to be blank */
  *whc    = 0.0;

  temp_whc  = 0.0;

  *return_code = CONIFERS_SUCCESS;

  if( dfsi < 50.0 && ppsi < 47.0 )
    {
      if( *whc <= 0.0 )
	{
	  *whc = 2.7;
	}
      if( *precip <= 0.0 )
	{
	  *precip = 35.0;
	}
      return;
    }

  else if( dfsi >= 50.0 && dfsi <= 160)
    {
      temp_whc = 0.22 * dfsi - 8.2;      
               
      if( *whc <= 0.0 ) 
	{
	  *whc  = temp_whc;
	}
      if( *precip <= 0.0 )
	{
	  *precip = 35.0;
	}
      return;
    }
  else if (ppsi >= 47.0 && ppsi <= 151)  
    {
      temp_whc = 0.22 * ppsi * 1.06 - 8.20;
      if ( *whc <= 0.0)
	{
	  *whc = temp_whc;
	}
      if ( *precip <= 0.0 )
	{
	  *precip= 35.0;
	}
      return;
    }
  else
    {
      temp_whc = 27.0;
      if ( *whc <= 0.0)
	{
	  *whc= temp_whc;
	}
      if( * precip <= 0.0)
	{
	  *precip = 35.0;
	}
      return;
    }
  *return_code = FILL_SAMPLE_FAILED;
  return;
  /*  if we get to here we have problems  */   

}


void calc_plot_stats_2( 
		       unsigned long           *return_code,
		       unsigned long           n_species,
		       struct SPECIES_RECORD   *species_ptr,
		       unsigned long           n_coeffs,
		       struct COEFFS_RECORD    *coeffs_ptr,
		       unsigned long           n_plants,
		       struct PLANT_RECORD     *plants_ptr,
		       unsigned long           n_points,
		       struct PLOT_RECORD      *plots_ptr )
{
        
  unsigned long   i;
  unsigned long   j;


  unsigned long   first_idx;
  unsigned long   last_idx; 
  struct  PLOT_RECORD     *plot_ptr;
  struct  PLANT_RECORD    *p_ptr; 
  unsigned long   n_plants_found;
  struct  COEFFS_RECORD   *c_ptr;

  /* go through the tree array and only tally the basal area  */
  /* and expf for those plants that are not shrubs            */
  *return_code = CONIFERS_ERROR;

  /* sort the tree list just to be safe   */
  /* take the performance hit...          */
  qsort(   plants_ptr, 
	   n_plants, 
	   sizeof( struct PLANT_RECORD ), 
	   compare_trees_by_plot ); 

  /* iterate through the plot and null out the values         */
  /* that will be calculated in the function, which should be */
  /* all of them                                              */
  plot_ptr = &plots_ptr[0];
  for( i = 0; i < n_points; i++, plot_ptr++ )
    {
      /* these are temp variables */
      plot_ptr->shrub_pct_cover   = 0.0;    /*  crown ratio calc                */
      plot_ptr->shrub_mean_height = 0.0;    /*  ditto                           */
      plot_ptr->shrub_expf        = 0.0;    /*  stems per acre for shrubs       */
      plot_ptr->basal_area        = 0.0;    /*  temp variable                   */
      plot_ptr->d6_area           = 0.0;    /*  temp variable                   */
      plot_ptr->expf              = 0.0;    /*  basal area of largest stem      */
      plot_ptr->bh_expf           = 0.0;    /*  expf in hw+con trees > 4.5 ft   */
      plot_ptr->qmd               = 0.0;    /*  quadratic mean diameter         */
      plot_ptr->sdi               = 0.0;    /*  stand density index             */
      plot_ptr->crown_area        = 0.0;    /*  crown area                      */
      plot_ptr->ccf               = 0.0;    /*  crown competition factor        */                                    
      plot_ptr->ba_c              = 0.0;    /*  b.h. basal area in conifers     */
      plot_ptr->ba_h              = 0.0;    /*  b.h. basal area in hardwoods    */
      plot_ptr->d6ba_c            = 0.0;    /*  basal ba in conifers            */
      plot_ptr->d6ba_h            = 0.0;    /*  basal ba in hardwoods           */
      plot_ptr->d6ba_s            = 0.0;    /*  basal area in shrubs            */                                    
      plot_ptr->ca_c              = 0.0;    /*  crown area in conifers          */
      plot_ptr->ca_h              = 0.0;    /*  crown area in hardwoods         */
      plot_ptr->ca_s              = 0.0;    /*  crown area in shrubs            */
      plot_ptr->hann_wang_x0      = 0.0;    /*  initial trajectory hann & wang  */

      /* added for CONIFERS_CIPS */
	plot_ptr->d12ba_c      = 0.0;
	  

	  /* clean out the old values */
	  memset( plot_ptr->bait, 0, sizeof( double ) * PLANT_TYPES * AIT_SIZE );
	  memset( plot_ptr->cait, 0, sizeof( double ) * PLANT_TYPES * AIT_SIZE );
      memset( plot_ptr->bal,  0, sizeof( double ) * PLANT_TYPES * AIT_SIZE );
    
	  /* the plants have to be sorted by plot */
	  get_plant_indecies_for_plot(  return_code,
					plot_ptr,
					n_plants,
					plants_ptr,
					&first_idx,
					&last_idx,
					&n_plants_found );
	  /* need_error_trap_here */

	  //p_ptr = &plants_ptr[first_idx+1];
	  p_ptr = &plants_ptr[first_idx];
	  for( j = first_idx; j <= last_idx; j++, p_ptr++ )
	    {
	      c_ptr = &coeffs_ptr[species_ptr[p_ptr->sp_idx].fsp_idx];

	      /* this is where the set_in_taller_attribs function is called */
        
	      set_in_taller_attribs( p_ptr, c_ptr, plot_ptr );
          /* added for CIPS variant August 1, 2014 by MWR*/
          set_in_larger_attribs( p_ptr, c_ptr, plot_ptr );

	      /* take care of the shrub information */
	      switch( c_ptr->type )
		{
		case CONIFER:
		  plot_ptr->ba_c         += p_ptr->basal_area * p_ptr->expf;
		  plot_ptr->ca_c         += p_ptr->crown_area * p_ptr->expf;
		  plot_ptr->d6ba_c       += p_ptr->d6_area    * p_ptr->expf;

          /* added for CONIFERS_CIPS */
          plot_ptr->d12ba_c    += p_ptr->d12_area * p_ptr->expf;

		  if(p_ptr->tht > 4.5 )
		    {
		      plot_ptr->bh_expf  += p_ptr->expf;
		    }


		  plot_ptr->basal_area   += p_ptr->basal_area * p_ptr->expf;
		  plot_ptr->ccf          += p_ptr->max_crown_width *
		    p_ptr->max_crown_width * 
		    p_ptr->expf;



          break;

		case HARDWOOD:
		  plot_ptr->ba_h         += p_ptr->basal_area * p_ptr->expf;
		  plot_ptr->ca_h         += p_ptr->crown_area * p_ptr->expf;
		  plot_ptr->d6ba_h       += p_ptr->d6_area    * p_ptr->expf;

		  if(p_ptr->tht > 4.5 )
		    {
		      plot_ptr->bh_expf  += p_ptr->expf;
		    }

		  plot_ptr->basal_area   += p_ptr->basal_area * p_ptr->expf;
		  plot_ptr->ccf          += p_ptr->max_crown_width *
		    p_ptr->max_crown_width * 
		    p_ptr->expf;
		  break;

		case SHRUB:
		  plot_ptr->shrub_pct_cover   += p_ptr->crown_area * p_ptr->expf;
		  plot_ptr->shrub_mean_height += p_ptr->tht        * p_ptr->expf;
		  plot_ptr->ca_s              += p_ptr->crown_area * p_ptr->expf;
		  plot_ptr->d6ba_s            += p_ptr->d6_area    * p_ptr->expf;
		  plot_ptr->shrub_expf        += p_ptr->expf;
		  break;

		  /* the forbs */
		case FORB:
		  /* and do what here ? */
		  break;
		}

	      plot_ptr->expf              += p_ptr->expf;
	      plot_ptr->d6_area           += p_ptr->d6_area * p_ptr->expf;
	      plot_ptr->crown_area        += p_ptr->crown_area * p_ptr->expf;
		
		  }



	  if(plot_ptr->shrub_expf>0.0)
	    {
	      plot_ptr->shrub_pct_cover   *= (100.0 / SQ_FT_PER_ACRE);
	      plot_ptr->shrub_mean_height /= plot_ptr->shrub_expf;
	    }

	  plot_ptr->ccf               *= CCF_CONST_I;   

	  if(plot_ptr->bh_expf > 0.0)
	    {
	      plot_ptr->qmd = 
		sqrt( (plot_ptr->basal_area / plot_ptr->bh_expf) / FC_I );
	      plot_ptr->sdi = 
		plot_ptr->bh_expf * pow( plot_ptr->qmd * 0.1, REINEKE_B1 );
	    }
	  else
	    {
	      plot_ptr->qmd = 0.0;
	      plot_ptr->sdi = 0.0;
	    }

	  if( plot_ptr->shrub_pct_cover < 0.0 )
	    {
	      plot_ptr->shrub_pct_cover = 0.0;
	    }

	  if( plot_ptr->shrub_mean_height < 0.0 )
	    {
	      plot_ptr->shrub_mean_height = 0.0;
	    }
	  //continue;

    }

  *return_code = CONIFERS_SUCCESS;
}



static int compare_trees_by_plot(
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
      return 0;
    }


}



static int compare_summaries_by_code( 
				     const void *ptr1, 
				     const void *ptr2 )
{
  struct SUMMARY_RECORD   *pt1_ptr;
  struct SUMMARY_RECORD   *pt2_ptr;

  pt1_ptr = (struct SUMMARY_RECORD*)ptr1;
  pt2_ptr = (struct SUMMARY_RECORD*)ptr2;
    
  //return strcmp( pt1_ptr->code, pt2_ptr->code );

  if( pt1_ptr->code > pt2_ptr->code )
    {
      return 1;
    }
  if( pt1_ptr->code < pt2_ptr->code )
    {
      return -1;
    }
  else
    {
      return 0;
    }

}




/********************************************************************************/
/*  Functional species summary functions                                        */
/********************************************************************************/

/* this allocates an array of summary records   */
/*  and fills in the species codes for summary */
struct SUMMARY_RECORD *build_fsp_summaries( 
					   unsigned long           *return_code,
					   unsigned long           n_species,
					   struct SPECIES_RECORD   *species_ptr,
					   unsigned long           n_coeffs,
					   struct COEFFS_RECORD    *coeffs_ptr,
					   unsigned long           n_plants,
					   struct PLANT_RECORD     *plants_ptr,
					   unsigned long           *n_fsp_in_sample )
{
    
  unsigned long           i;
  struct  PLANT_RECORD    *plant_ptr;
  //struct  COEFFS_RECORD   *fsp_ptr;
  struct  COEFFS_RECORD   *c_ptr;
  struct  SUMMARY_RECORD  *fsp_sum_ptr;
  struct  SUMMARY_RECORD  *temp_sum_ptr;

  fsp_sum_ptr         = NULL;
  *n_fsp_in_sample    = 0;

  /* start filling in the array from the plant list               */
  /* first find the number of unique species in the plant list    */
  plant_ptr = &plants_ptr[0];
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        
      c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

      /* see if the speices for the current plant is already      */
      /* in the sum_ptr array                                     */
      /* if it is, continue, otherwise add it                     */
      /* man is this a bad function name                          */
      temp_sum_ptr = get_summary_from_code(   *n_fsp_in_sample,
					      fsp_sum_ptr,
					      c_ptr->idx );

      /* if it's found, then the entry is already in the sum_ptr  */
      /* if the temp_sum_ptr == NULL, then you coulnd't find it   */
      /* so add it and increment the number entered               */
      if( temp_sum_ptr != NULL )
	{
	  continue;
	}
      else
	{
	  /* if it's not in the sample, you'll need to add a new  */
	  /* element by realloc-ing the array and adding the new  */
	  /* value to the end, and resorting                      */
	  fsp_sum_ptr = 
	    (struct SUMMARY_RECORD *)realloc(   
					     fsp_sum_ptr,
					     *n_fsp_in_sample * sizeof( struct SUMMARY_RECORD ) + 
					     sizeof( struct SUMMARY_RECORD ) );

	  fsp_sum_ptr[(*n_fsp_in_sample)++].code = c_ptr->idx;
            
	  /* now sort so you can use bsearch() to see if the entry    */
	  /* has already been added                                   */
	  qsort(  fsp_sum_ptr, 
		  *n_fsp_in_sample, 
		  sizeof( struct SUMMARY_RECORD ), 
		  compare_summaries_by_code );  
	}

    }

  /* now sort the fsp summaries by fsp */
  qsort(  fsp_sum_ptr, 
	  *n_fsp_in_sample, 
	  sizeof( struct SUMMARY_RECORD ), 
	  compare_summaries_by_code );  

  /* return the sumamry array pointer             */
  /* set the return_code for the calling function */
  /* and resort the tree list                     */
  *return_code = CONIFERS_SUCCESS;


  return fsp_sum_ptr;

}



/* this function updates the funstional species sumamries */
void update_fsp_summaries( 
			  unsigned long           *return_code,
			  unsigned long           n_species,
			  struct SPECIES_RECORD   *species_ptr,
			  unsigned long           n_coeffs,
			  struct COEFFS_RECORD    *coeffs_ptr,
			  unsigned long           n_plants,
			  struct PLANT_RECORD     *plants_ptr,
			  unsigned long           n_points,
			  unsigned long           n_summaries,
			  struct SUMMARY_RECORD   *summaries_ptr )
{
    
  unsigned long           i;
  struct  COEFFS_RECORD   *c_ptr;
  struct  PLANT_RECORD    *plant_ptr;
  struct  SUMMARY_RECORD  *sum_ptr;
  unsigned long   temp_code;
  double  temp_volume;                    /*  MOD025  */
  double  temp_biomass;                 /*  MOD026  */  
  double  hdr;                            /*  MOD042  */

  /* iterate through the sumamries and null out the values    */
  /* that will be calculated in the function, which should be */
  /* all of them                                              */
  sum_ptr = &summaries_ptr[0];
  for( i = 0; i < n_summaries; i++, sum_ptr++ )
    {
      /* store the sp_code in a temp and memset the struct    */
      /* copy the sp_code back in so you don't have to        */
      /* keep nulling out the struct by hand                  */
      temp_code = sum_ptr->code;
      memset( sum_ptr, 0, sizeof( struct  SUMMARY_RECORD ) );
      sum_ptr->code = temp_code;
    }

  /* go through the loop AGAIN! to get the min/max values for */
  /* min_dbh/max_dbh and height                               */
  plant_ptr = &plants_ptr[0];
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {

      c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

      /* now get the summary code for the plant record */            
      sum_ptr = get_summary_from_code(    n_summaries,
					  summaries_ptr,
					  c_ptr->idx );
      if( sum_ptr == NULL )
	{
	  *return_code = INVALID_SP_CODE;
	  return;
	}
      else
	{
	  /* calc the min dbh values */
	  /* start the array by setting the values to the max values */
	  if( plant_ptr->dbh > sum_ptr->min_dbh && plant_ptr->dbh > 0.0 )
	    {
	      sum_ptr->min_dbh = plant_ptr->dbh;
	    }

	  /* calc the min height values */
	  if( plant_ptr->tht > sum_ptr->min_height && plant_ptr->tht > 0.0 )
	    {
	      sum_ptr->min_height = plant_ptr->tht;
	    }

	  /* MOD042 */
	  /* calc the min h/d values */
	  hdr = plant_ptr->tht / plant_ptr->d6;
	  if( hdr > sum_ptr->min_hd6_ratio && hdr > 0.0 )
	    {
	      sum_ptr->min_hd6_ratio = hdr;
	    }

	  //continue;
	}
    }


  /* start filling in the array from the plant list               */
  /* first find the number of unique species in the plant list    */
  plant_ptr = &plants_ptr[0];
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {

      c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];
        
      /* see if the speices for the current plant is already      */
      /* in the sum_ptr array                                     */
      /* if it is, continue, otherwise add it                     */
      sum_ptr = get_summary_from_code(    n_summaries,
					  summaries_ptr,
					  c_ptr->idx );
         
      /* if it's found, then the entry is already in the sum_ptr  */
      /* if the temp_sum_ptr == NULL, then you coulnd't find it   */
      /* and something's very wrong, otherwise you're cool        */
      if( sum_ptr == NULL )
	{
	  *return_code = INVALID_SP_CODE;
	  return;
	}
      else
	{
	  sum_ptr->expf           += plant_ptr->expf;

	  /* add the ccf values for the sample */
	  c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

	  /*  MOD034 */
	  if( is_tree( c_ptr ) )
	    {         
	      sum_ptr->ccf            += CCF_CONST_I *
									 plant_ptr->max_crown_width * 
									 plant_ptr->max_crown_width *
									 plant_ptr->expf;
	      sum_ptr->cr             += plant_ptr->cr * plant_ptr->expf;
	      sum_ptr->tree_expf      += plant_ptr->expf;

	      hdr = plant_ptr->tht / plant_ptr->d6;
	      sum_ptr->mean_hd6_ratio  += hdr * plant_ptr->expf;

	      if( hdr > sum_ptr->min_hd6_ratio )
		  {
		    sum_ptr->max_hd6_ratio = hdr;
		  }

	      if( hdr < sum_ptr->min_hd6_ratio )
		  {
		    sum_ptr->min_hd6_ratio = hdr;
		  }

	      if(c_ptr->type == CONIFER)
		  {
		    sum_ptr->con_tpa += plant_ptr->expf;
		  }
	    }


	  /* MOD012 */
	  /*  calculate the expf for trees        */
	  /*  above bh (conifs and hwoods only)   */
	  if ( plant_ptr->tht > 4.5 )
	    {           
	      /* only sum up the values for the trees that    */
	      /* are over 4.5 feet  tall                      */
	      /*  MOD034  */
	      if( is_tree( c_ptr ) )
		{         
		  sum_ptr->bh_expf    += plant_ptr->expf;
		  sum_ptr->basal_area += plant_ptr->basal_area * plant_ptr->expf;
		  calc_volume(return_code,
			      plant_ptr->tht,
			      plant_ptr->dbh,
			      &temp_volume,
			      c_ptr->cfvolume4);
		  /* need_error_trap_here */
		  sum_ptr->cfvolume4  += temp_volume * plant_ptr->expf;
		}
	    }
            
	  calc_biomass(return_code, /* MOD026 call biomass function */ 
		       plant_ptr->tht,
		       plant_ptr->d6,
		       plant_ptr->crown_width,
		       plant_ptr->dbh,
		       &temp_biomass,
		       c_ptr->biomass);
	  /* need_error_trap_here */
	  sum_ptr->biomass      += temp_biomass * plant_ptr->expf; 
	  sum_ptr->mean_height    += plant_ptr->tht * plant_ptr->expf;
	  sum_ptr->crown_area     += plant_ptr->crown_area * plant_ptr->expf;
	  sum_ptr->pct_cover      += plant_ptr->crown_area * plant_ptr->expf;

	  /* calc the min/max height values */
	  if( plant_ptr->tht < sum_ptr->min_height && plant_ptr->tht > 0.0 )
	    {
	      sum_ptr->min_height = plant_ptr->tht;
	    }
	  if( plant_ptr->tht > sum_ptr->max_height )
	    {
	      sum_ptr->max_height = plant_ptr->tht;
	    }


	  /* calc the min/max dbh values */
	  if( plant_ptr->dbh < sum_ptr->min_dbh && 
	      plant_ptr->dbh > 0.0 )
	    {
	      sum_ptr->min_dbh = plant_ptr->dbh;
	    }
	  if( plant_ptr->dbh > sum_ptr->max_dbh )
	    {
	      sum_ptr->max_dbh = plant_ptr->dbh;
	    }

	}

    }

  /* adjust the species summaries by the number of plots */
  sum_ptr = &summaries_ptr[0];
  for( i = 0; i < n_summaries; i++, sum_ptr++ )
    {
      sum_ptr->cr             /= sum_ptr->tree_expf;  /* MOD034 */
        
      sum_ptr->mean_hd6_ratio /= sum_ptr->expf;   /* MOD042 */
      sum_ptr->tree_expf      /= n_points;
      sum_ptr->mean_height    /= sum_ptr->expf;        
      sum_ptr->bh_expf        /= n_points;
      sum_ptr->basal_area     /= n_points;
      sum_ptr->expf           /= n_points;
      sum_ptr->crown_area     /= n_points;
      sum_ptr->ccf            /= n_points;
      sum_ptr->cfvolume4      /= n_points;  /* MOD025 */
      sum_ptr->biomass      /= n_points;    /* MOD026 */
      sum_ptr->pct_cover      /= (n_points * SQ_FT_PER_ACRE);
      sum_ptr->pct_cover      *= 100.0;
      sum_ptr->con_tpa        /= n_points;

      /* MOD042 */
      if( sum_ptr->expf <= 0.0 )
	{
	  sum_ptr->min_hd6_ratio  = 0.0;   /* MOD042 */
	  sum_ptr->mean_hd6_ratio = 0.0;   /* MOD042 */
	  sum_ptr->max_hd6_ratio  = 0.0;   /* MOD042 */
	}
        
      /* MOD011 */
      if( sum_ptr->bh_expf > 0.0 )
	{
	  sum_ptr->qmd = sqrt( (sum_ptr->basal_area / sum_ptr->bh_expf) / FC_I );
	  sum_ptr->sdi = sum_ptr->bh_expf * pow( sum_ptr->qmd * 0.1, REINEKE_B1 );
	}
      else
	{
	  sum_ptr->qmd = 0.0;
	  sum_ptr->sdi = 0.0;
	}            

    }

  *return_code = CONIFERS_SUCCESS;
}

/********************************************************************************/
/*  Species summary functions                                                   */
/********************************************************************************/

/* this allocates an array of summary records   */
/*  and fills in the species codes for summary */
struct SUMMARY_RECORD *build_species_summaries( 
					       unsigned long           *return_code,
					       unsigned long           n_species,
					       struct SPECIES_RECORD   *species_ptr,
					       unsigned long           n_plants,
					       struct PLANT_RECORD     *plants_ptr,
					       unsigned long           *n_sp_in_sample )
{
   
  unsigned long           i;
  struct  PLANT_RECORD    *plant_ptr = NULL;
  struct  SPECIES_RECORD  *sp_ptr = NULL;
  struct  SUMMARY_RECORD  *sp_sum_ptr = NULL;
  struct  SUMMARY_RECORD  *temp_sum_ptr = NULL;

  *n_sp_in_sample    = 0;
  *return_code = CONIFERS_SUCCESS;

  /* start filling in the array from the plant list               */
  /* first find the number of unique species in the plant list    */
  plant_ptr = &plants_ptr[0];
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        
      /* get the species index for the plant */
      sp_ptr = &species_ptr[plant_ptr->sp_idx];

      /* see if the speices for the current plant is already      */
      /* in the sum_ptr array                                     */
      /* if it is, continue, otherwise add it                     */
      /* man is this a bad function name                          */
      temp_sum_ptr = get_summary_from_code( *n_sp_in_sample,
					    sp_sum_ptr,
					    sp_ptr->idx );
      /* if it's found, then the entry is already in the sum_ptr  */
      /* if the temp_sum_ptr == NULL, then you coulnd't find it   */
      /* so add it and increment the number entered               */
      if( temp_sum_ptr != NULL )
	{
	  continue;
	}
      else
	{
	  /* if it's not in the sample, you'll need to add a new  */
	  /* element by realloc-ing the array and adding the new  */
	  /* value to the end, and resorting                      */
	  sp_sum_ptr = 
	    (struct SUMMARY_RECORD *)realloc(   
					     sp_sum_ptr,
					     *n_sp_in_sample * sizeof( struct SUMMARY_RECORD ) + 
					     sizeof( struct SUMMARY_RECORD ) );

	  /* MOD028 */
	  /* NULL out the entry */
	  memset( &sp_sum_ptr[(*n_sp_in_sample)], 
		  0, 
		  sizeof( struct SUMMARY_RECORD ) );

	  //strcpy( sp_sum_ptr[(*n_sp_in_sample)++].code, sp_ptr->sp_code );
	  sp_sum_ptr[(*n_sp_in_sample)++].code = sp_ptr->idx;
            
	  /* now sort so you can use bsearch() to see if the entry    */
	  /* has already been added                                   */
	  qsort(  sp_sum_ptr, 
		  *n_sp_in_sample, 
		  sizeof( struct SUMMARY_RECORD ), 
		  compare_summaries_by_code );
	}

    }



  /* now sort the fsp summaries by fsp */
  qsort(  sp_sum_ptr, 
	  *n_sp_in_sample, 
	  sizeof( struct SUMMARY_RECORD ), 
	  compare_summaries_by_code );  

  /* return the sumamry array pointer             */
  /* set the return_code for the calling function */
  /* and resort the tree list                     */
  *return_code = CONIFERS_SUCCESS;

  return sp_sum_ptr;

}


void update_species_summaries( 
			      unsigned long           *return_code,
			      unsigned long           n_species,
			      struct SPECIES_RECORD   *species_ptr,
			      unsigned long           n_coeffs,
			      struct COEFFS_RECORD    *coeffs_ptr,
			      unsigned long           n_plants,
			      struct PLANT_RECORD     *plants_ptr,
			      unsigned long           n_points,
			      unsigned long           n_summaries,
			      struct SUMMARY_RECORD   *summaries_ptr )
{
    
  unsigned long            i;
  unsigned long            j;
  struct  COEFFS_RECORD    *c_ptr;
  struct  PLANT_RECORD     *plant_ptr;
  struct  SUMMARY_RECORD   *sum_ptr;
  unsigned long            temp_code;
  double                   temp_volume;                
  double                   temp_biomass;               
  double                   hdr;                        
  struct                   HTN_RECORD *ht_n;
  double                   tht40;
  double                   sum_exp;


  /*note: in this routine ht_n->is_tree is assigned the code from the summaries_ptr if its a conifer */
  /* elsewhere it actually holds an indicator of is.tree either way it is an integer....             */
  ht_n = (struct HTN_RECORD*)calloc( n_plants, sizeof( struct HTN_RECORD ) );

  if( ht_n == NULL )
    {
      *return_code    = FAILED_MEMORY_ALLOC;
      return;
    }

  /* iterate through the sumamries and null out the values    */
  /* that will be calculated in the function, which should be */
  /* all of them                                              */
  sum_ptr = &summaries_ptr[0];
  for( i = 0; i < n_summaries; i++, sum_ptr++ )
    {
      /* store the sp_code in a temp and memset the struct    */
      /* copy the sp_code back in so you don't have to        */
      /* keep nulling out the struct by hand                  */
      temp_code = sum_ptr->code;
      memset( sum_ptr, 0, sizeof( struct  SUMMARY_RECORD ) );
      sum_ptr->code = temp_code;
    }

  /* go through the loop AGAIN! to get the min/max values for */
  /* min_dbh/max_dbh and height                               */
  plant_ptr = &plants_ptr[0];
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
      c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];


      sum_ptr = get_summary_from_code(    n_summaries,
					  summaries_ptr,
					  //plant_ptr->sp_code );
					  plant_ptr->sp_idx );
      if( sum_ptr == NULL )
	{
	  *return_code = INVALID_SP_CODE;
	  free(ht_n);
	  return;
	}
      else
	{
	  if (c_ptr->type == CONIFER) 
	    {
	      ht_n[i].is_tree = (int)sum_ptr->code;
	      ht_n[i].tht = plant_ptr->tht;
	      ht_n[i].expf = plant_ptr->expf/n_points;
	    }

	  if( plant_ptr->dbh > sum_ptr->min_dbh && plant_ptr->dbh > 0.0 )
	    {
	      sum_ptr->min_dbh = plant_ptr->dbh;
	    }

	  /* calc the min height values */
	  if( plant_ptr->tht > sum_ptr->min_height && plant_ptr->tht > 0.0 )
	    {
	      sum_ptr->min_height = plant_ptr->tht;
	    }

	  /* MOD042 */
	  /* calc the min h/d values */
	  if( is_tree( c_ptr ) && plant_ptr->d6 > 0.0 )
	    {
	      hdr = plant_ptr->tht / plant_ptr->d6;

	      if( hdr > sum_ptr->min_hd6_ratio && hdr > 0.0 )
		{
		  sum_ptr->min_hd6_ratio = hdr;
		}
	    }
	}
    }

  /* start filling in the array from the plant list               */
  /* first find the number of unique species in the plant list    */
  plant_ptr = &plants_ptr[0];
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
      /* see if the speices for the current plant is already      */
      /* in the sum_ptr array                                     */
      /* if it is, continue, otherwise add it                     */
      sum_ptr = get_summary_from_code(    n_summaries,
					  summaries_ptr,
					  //plant_ptr->sp_code );
					  plant_ptr->sp_idx );
         
      /* if it's found, then the entry is already in the sum_ptr  */
      /* if the temp_sum_ptr == NULL, then you coulnd't find it   */
      /* and something's very wrong, otherwise you're cool        */
      if( sum_ptr == NULL )
	{
	  *return_code = INVALID_SP_CODE;
	  free(ht_n);
	  return;
	}
      else
	{
	  sum_ptr->expf           += plant_ptr->expf;
              
	  c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

	  /*   MOD034   */
	  if( is_tree( c_ptr ) )
	    {         
	      sum_ptr->ccf            +=  CCF_CONST_I *
		plant_ptr->max_crown_width * 
		plant_ptr->max_crown_width *
		plant_ptr->expf;
	      sum_ptr->cr             += plant_ptr->cr * plant_ptr->expf;
	      sum_ptr->tree_expf      += plant_ptr->expf;

	      hdr = plant_ptr->tht / plant_ptr->d6;
	      sum_ptr->mean_hd6_ratio  += hdr * plant_ptr->expf;

	      if( hdr > sum_ptr->max_hd6_ratio )
		{
		  sum_ptr->max_hd6_ratio = hdr;
		}

	      if( hdr < sum_ptr->min_hd6_ratio )
		{
		  sum_ptr->min_hd6_ratio = hdr;
		}
	      if(c_ptr->type == CONIFER)
		{
		  sum_ptr->con_tpa += plant_ptr->expf;
		}
	    }


	  /* MOD012 */
	  /*  calculate the expf for trees        */
	  /*  above bh (conifs and hwoods only)   */
	  if ( plant_ptr->tht > 4.5 )
	    {           
	      /* only sum up the values for the trees that    */
	      /* are over 4.5 feet  tall                      */
	      if( is_tree( c_ptr ) )
		{         
		  sum_ptr->bh_expf    += plant_ptr->expf;
		  sum_ptr->basal_area += plant_ptr->basal_area * plant_ptr->expf;
                    
		  /* MOD025 */
		  calc_volume(return_code,  
			      plant_ptr->tht,
			      plant_ptr->dbh,
			      &temp_volume,
			      c_ptr->cfvolume4);

		    /* need_error_trap_here */
		    sum_ptr->cfvolume4  += temp_volume * plant_ptr->expf;
		}
	    }

	  calc_biomass(return_code, /* MOD026 call biomass function */ 
		       plant_ptr->tht,
		       plant_ptr->d6,
		       plant_ptr->crown_width,
		       plant_ptr->dbh,
		       &temp_biomass,
		       c_ptr->biomass);
	  /* need_error_trap_here */
	  sum_ptr->biomass      += temp_biomass * plant_ptr->expf; 
	  sum_ptr->mean_height    += plant_ptr->tht * plant_ptr->expf;
	  sum_ptr->crown_area     += plant_ptr->crown_area * plant_ptr->expf;
	  sum_ptr->pct_cover      += plant_ptr->crown_area * plant_ptr->expf;

	  /* calc the min/max height values */
	  if( plant_ptr->tht < sum_ptr->min_height && plant_ptr->tht > 0.0 )
	    {
	      sum_ptr->min_height = plant_ptr->tht;
	    }
	  if( plant_ptr->tht > sum_ptr->max_height )
	    {
	      sum_ptr->max_height = plant_ptr->tht;
	    }

	  /* calc the min/max dbh values */
	  if( plant_ptr->dbh < sum_ptr->min_dbh && 
	      plant_ptr->dbh > 0.0 )
	    {
	      sum_ptr->min_dbh = plant_ptr->dbh;
	    }
	  if( plant_ptr->dbh > sum_ptr->max_dbh )
	    {
	      sum_ptr->max_dbh = plant_ptr->dbh;
	    }

	}
    }            /* end of second tree loop*/



  /* adjust the species summaries by the number of plots */
  sum_ptr = &summaries_ptr[0];
  for( i = 0; i < n_summaries; i++, sum_ptr++ )
    {
      sum_ptr->cr             /= sum_ptr->tree_expf;  /*  MOD034  */
      sum_ptr->mean_hd6_ratio /= sum_ptr->expf;       /*  MOD042  */
      sum_ptr->tree_expf      /= n_points;
      sum_ptr->mean_height    /= sum_ptr->expf;        
      sum_ptr->bh_expf        /= n_points;
      sum_ptr->basal_area     /= n_points;
      sum_ptr->expf           /= n_points;
      sum_ptr->crown_area     /= n_points;
      sum_ptr->ccf            /= n_points;
      sum_ptr->cfvolume4      /= n_points;            /*   MOD025   */
      sum_ptr->biomass        /= n_points;            /*   MOD026   */
      sum_ptr->pct_cover      /= (n_points * SQ_FT_PER_ACRE);
      sum_ptr->pct_cover      *= 100.0;
      sum_ptr->con_tpa        /= n_points;
                
      /* MOD042 */
      if( sum_ptr->expf <= 0.0 )
	{
	  sum_ptr->min_hd6_ratio  = 0.0;   /* MOD042 */
	  sum_ptr->mean_hd6_ratio = 0.0;   /* MOD042 */
	  sum_ptr->max_hd6_ratio  = 0.0;   /* MOD042 */
	}


      /* MOD011 */
      if( sum_ptr->bh_expf > 0.0 )
	{
	  sum_ptr->qmd = sqrt( (sum_ptr->basal_area / sum_ptr->bh_expf) / FC_I );
	  sum_ptr->sdi = sum_ptr->bh_expf * pow( sum_ptr->qmd * 0.1, REINEKE_B1 );
	}
      else
	{
	  sum_ptr->qmd = 0.0;
	  sum_ptr->sdi = 0.0;
	}            

      /* MOD042 */
      /* check for extreme values */
      if( sum_ptr->mean_height < 0.0 ) 
	{
	  sum_ptr->mean_height = 0.0;
	}
            
      if( sum_ptr->cr < 0.0 ) 
	{
	  sum_ptr->cr = 0.0;
	}

    }

  qsort(   ht_n, 
	   n_plants, 
	   sizeof( struct HTN_RECORD ), 
	   compare_htn_by_plant_tht_expf ); 

  tht40   = 0.0;   /* initialize ht40 */
  sum_exp = 0.0;   /* init sum exps for ht40*/

  sum_ptr = &summaries_ptr[0];
  for( j = 0; j < n_summaries; j++, sum_ptr++ )
    {
      tht40   = 0.0;   /* initialize ht40 */
      sum_exp = 0.0;   /* init sum exps for ht40*/
      for (i = 0; i < n_plants; i++)
	{
	  if(sum_exp < 40.0 && ht_n[i].is_tree == (int)sum_ptr->code)
	    {
	      if( sum_exp + ht_n[i].expf <= 40.0   )
		{
		  tht40   = tht40   + ht_n[i].tht * ht_n[i].expf;
		  sum_exp = sum_exp + ht_n[i].expf;
		}
	      else 
		{
		  tht40 = tht40   + ht_n[i].tht * (40 - sum_exp );
		  sum_exp = 40.0;
		}
	    }
	}
      if(sum_exp >0.0)
	{
	  sum_ptr->height_40=tht40/sum_exp;
	}
    }

  free(ht_n);

  *return_code = CONIFERS_SUCCESS;
}


static int compare_htn_by_plant_tht_expf( 
					 const void *ptr1, 
					 const void *ptr2 )
{
  struct HTN_RECORD   *pt1_ptr;
  struct HTN_RECORD   *pt2_ptr;

  pt1_ptr = (struct HTN_RECORD*)ptr1;
  pt2_ptr = (struct HTN_RECORD*)ptr2;

  if( pt1_ptr->is_tree < pt2_ptr->is_tree )
    {
      return 1;
    }
  if( pt1_ptr->is_tree > pt2_ptr->is_tree )
    {
      return -1;
    }
  else
    {
      if( pt1_ptr->tht < pt2_ptr->tht )
	{
	  return 1;
	}
      if( pt1_ptr->tht > pt2_ptr->tht )
	{
	  return -1;
	}
      else
	{
	  if( pt1_ptr->expf < pt2_ptr->expf )
	    {
	      return 1;
	    }
	  if( pt1_ptr->expf > pt2_ptr->expf )
	    {
	      return -1;
	    }
	  else
	    {
	      return 0;
	    }
	}
    }
}



/*  MOD025  */
/*  MOD003  */
/* this function updates a single struct SUMMARY_RECORD     */
/* that represents the summary record for all the species   */
/* in the sample. Note not all values will add up from the  */
/* component species found in the species_ptr array         */
void update_total_summaries( 
			    unsigned long           *return_code,
			    unsigned long           n_points,
			    unsigned long           n_plants,
			    unsigned long           n_species,
			    struct SPECIES_RECORD   *species_ptr,
			    unsigned long           n_coeffs,
			    struct COEFFS_RECORD    *coeffs_ptr,
			    struct PLANT_RECORD     *plants_ptr,
			    struct SUMMARY_RECORD   *sum_ptr )
{
    
  unsigned long           i;
  struct  PLANT_RECORD    *plant_ptr;
  struct  COEFFS_RECORD   *c_ptr;

  double                  max_sdi;
  double                    temp_volume;    /*  MOD025  */
  double                    temp_biomass;   /*  MOD026  */
  //   double                  *double_ptr;    /*  MOD043  */
  double                  hdr;            /*  MOD042  */
  double  tht40;
  double  sum_exp;
  double  temp_min_d;
  double  temp_min_h;
  struct HTN_RECORD *ht_n;



  /* iterate through the sumamries and null out the values    */
  /* that will be calculated in the function, which should be */
  /* all of them                                              */
  memset( sum_ptr, 0, sizeof( struct SUMMARY_RECORD ) );

  /* start filling in the array from the plant list               */
  /* first find the number of unique species in the plant list    */
  plant_ptr = &plants_ptr[0];

  // added this block of code April 2008 for initializing mins and maxs
  temp_min_d = 0.0;
  temp_min_h = 0.0;
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
      if((plant_ptr->tht > temp_min_h) && (plant_ptr->tht > 0.0))
	{
	  temp_min_h = plant_ptr->tht;
	}
      if((plant_ptr->dbh > temp_min_d) && (plant_ptr->dbh > 0.0))
	{
	  temp_min_d = plant_ptr->dbh;
	}
    }
  
  /* just start the min/max values on the first plant */
  sum_ptr->min_dbh = temp_min_d;
  sum_ptr->max_dbh = 0.0;

  sum_ptr->min_height = temp_min_h;
  sum_ptr->max_height = 0.0;

  /* just start the min/max values on the first plant */

  /* allocate an array of three doubles that's n_plants long */
  //*ht_n = (double*)calloc( n_plants, sizeof( double ) * 3 * n_plants );
  ht_n = (struct HTN_RECORD*)calloc( n_plants, sizeof( struct HTN_RECORD ) );

  if( ht_n == NULL )
    {
      *return_code    = FAILED_MEMORY_ALLOC;
      return;
    }

  plant_ptr = &plants_ptr[0];

  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
      sum_ptr->expf           += plant_ptr->expf;

      c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

      /* only sum up the values for the trees that    */
      /* are over 4.5 feet  tall                      */
      if( plant_ptr->tht > 4.5 )
	{
	  /*  MOD034  */
	  if( is_tree( c_ptr ) )
	    {         
	      sum_ptr->bh_expf    += plant_ptr->expf;
	      sum_ptr->basal_area += plant_ptr->basal_area * plant_ptr->expf;
                
	      /*  MOD025  */
	      calc_volume(    return_code, 
			      plant_ptr->tht,
			      plant_ptr->dbh,
			      &temp_volume,
			      c_ptr->cfvolume4 );
	      /* need_error_trap_here */

	      sum_ptr->cfvolume4  += temp_volume * plant_ptr->expf;
	    }
	}

      /* MOD042 */
      /* calc the min h/d values */
      if( is_tree( c_ptr ) && plant_ptr->d6 > 0.0 )
	{
	  hdr = plant_ptr->tht / plant_ptr->d6;

	  if( hdr > sum_ptr->min_hd6_ratio && hdr > 0.0 )
	    {
	      sum_ptr->min_hd6_ratio = hdr;
	    }
	}
      /* this is to calculate the height 40*/
      if (c_ptr->type == CONIFER) 
	{
	  ht_n[i].is_tree = is_tree( c_ptr );
	  ht_n[i].tht = plant_ptr->tht;
	  ht_n[i].expf = plant_ptr->expf/n_points;
	}

      calc_biomass(return_code, /* MOD026 call biomass function */ 
		   plant_ptr->tht,
		   plant_ptr->d6,
		   plant_ptr->crown_width,
		   plant_ptr->dbh,
		   &temp_biomass,
		   c_ptr->biomass);
      /* need_error_trap_here */
      sum_ptr->biomass          += temp_biomass * plant_ptr->expf; 

      sum_ptr->crown_area     += plant_ptr->crown_area * plant_ptr->expf;
      sum_ptr->pct_cover      += plant_ptr->crown_area * plant_ptr->expf;

      /* add the ccf values for the sample */
      if( is_tree( c_ptr ) )
	{
	  sum_ptr->mean_height    += plant_ptr->tht * plant_ptr->expf;
	  sum_ptr->ccf            +=  CCF_CONST_I *
	    plant_ptr->max_crown_width * 
	    plant_ptr->max_crown_width *
	    plant_ptr->expf;
	  sum_ptr->cr             += plant_ptr->cr * plant_ptr->expf;
	  sum_ptr->tree_expf      += plant_ptr->expf;
	  if(c_ptr->type == CONIFER)
	    {
	      sum_ptr->con_tpa += plant_ptr->expf;
	    }
	}

      /* calc the min/max height values */
      if( is_tree(c_ptr) )
	{
	  if( plant_ptr->tht < sum_ptr->min_height && plant_ptr->tht > 0.0 )
	    {
	      sum_ptr->min_height = plant_ptr->tht;
	    }
	  if( plant_ptr->tht > sum_ptr->max_height )
	    {
	      sum_ptr->max_height = plant_ptr->tht;
	    }
        
	  /* calc the min/max dbh values */
	  if( plant_ptr->dbh < sum_ptr->min_dbh && 
	      plant_ptr->dbh > 0.0 )
	    {
	      sum_ptr->min_dbh = plant_ptr->dbh;
	    }
	  if( plant_ptr->dbh > sum_ptr->max_dbh )
	    {
	      sum_ptr->max_dbh = plant_ptr->dbh;
	    }
	}
    }

  /*  MOD024  */
  /*  MOD034  */
  /* adjust the species summaries by the number of plots */
  sum_ptr->cr             /= sum_ptr->tree_expf;
  sum_ptr->mean_hd6_ratio /= sum_ptr->expf;   /* MOD042 */
  sum_ptr->mean_height    /= sum_ptr->tree_expf;
  sum_ptr->tree_expf      /= n_points;
  sum_ptr->bh_expf        /= n_points;
  sum_ptr->basal_area     /= n_points;
  sum_ptr->expf           /= n_points;
  sum_ptr->crown_area     /= n_points;
  sum_ptr->ccf            /= n_points;  /*  MOD035   */
  sum_ptr->cfvolume4      /= n_points;  /*  MOD025   */
  sum_ptr->biomass        /= n_points;  /*  MOD026   */
  sum_ptr->pct_cover      /= (n_points * SQ_FT_PER_ACRE);
  sum_ptr->pct_cover      *= 100.0;
  sum_ptr->con_tpa         /= n_points;

  /* MOD042 */
  if( sum_ptr->expf <= 0.0 )
    {
      sum_ptr->min_hd6_ratio  = 0.0;   /* MOD042 */
      sum_ptr->mean_hd6_ratio = 0.0;   /* MOD042 */
      sum_ptr->max_hd6_ratio  = 0.0;   /* MOD042 */
    }


  /* MOD011 */
  if( sum_ptr->bh_expf > 0.0 )
    {
      sum_ptr->qmd = sqrt( (sum_ptr->basal_area / sum_ptr->bh_expf) / FC_I );
      sum_ptr->sdi = sum_ptr->bh_expf * pow( sum_ptr->qmd * 0.1, REINEKE_B1 );
    }
  else
    {
      sum_ptr->qmd = 0.0;
      sum_ptr->sdi = 0.0;
    }            

  /* check for extream values */
  if( sum_ptr->mean_height < 0.0 ) 
    {
      sum_ptr->mean_height = 0.0;
    }

  if( sum_ptr->cr < 0.0 ) 
    {
      sum_ptr->cr = 0.0;
    }

  /* default values */
  sum_ptr->curtis_rd      = 0.0;
  sum_ptr->rel_density    = 0.0;
  max_sdi                 = 0.0;
  if( sum_ptr->qmd > 0.0  && sum_ptr->basal_area >0.0)
    {
      sum_ptr->curtis_rd = sum_ptr->basal_area / sqrt( sum_ptr->qmd );
      
      calc_max_sdi( return_code, 
		    n_species,
		    species_ptr,
		    n_coeffs,
		    coeffs_ptr,
		    n_plants,
		    plants_ptr,
		    n_points,
		    &max_sdi );
    }

  sum_ptr->sdimax = max_sdi;
  /* now calculate the REINEKE relative density */
  if( max_sdi > 0.0 )
    {
      sum_ptr->rel_density = sum_ptr->sdi / max_sdi;
    }
   
  /* sort for ht 40 */
  qsort(   ht_n, 
	   n_plants, 
	   sizeof( struct HTN_RECORD ), 
	   compare_htn_by_plant_tht_expf ); 

  tht40   = 0.0;   /* initialize ht40 */
  sum_exp = 0.0;   /* init sum exps for ht40*/

  for (i = 0; i < n_plants; i++)
    {
      /* first accumulate height of 40 largest conifers...*/
      if(sum_exp < 40.0 && ht_n[i].is_tree ==1)
	{
	  if( sum_exp + ht_n[i].expf <= 40.0   )
	    {
	      tht40   = tht40   + ht_n[i].tht * ht_n[i].expf;
	      sum_exp = sum_exp + ht_n[i].expf;
	    }
	  else 
	    {
	      tht40 = tht40   + ht_n[i].tht * (40 - sum_exp );
	      sum_exp = 40.0;
	    }
	}
    }

  free(ht_n);
  if(sum_exp > 0.0)
    {
      sum_ptr->height_40=tht40/sum_exp;
    }

  *return_code = CONIFERS_SUCCESS;

}





struct SUMMARY_RECORD *get_summary_from_code(
					     unsigned  long          n_summaries,
					     struct SUMMARY_RECORD   *summaries_ptr,
					     unsigned long           code )
{

  struct SUMMARY_RECORD *entry;
  struct SUMMARY_RECORD key;

  /* now use b_search to find the coeffs entry */
  memset( &key, 0, sizeof( struct SUMMARY_RECORD ) );
  key.code = code;

  /* find the matching summary record */
  entry = (struct SUMMARY_RECORD*)bsearch( 
					  &key, 
					  summaries_ptr, 
					  (size_t)n_summaries,
					  sizeof( struct SUMMARY_RECORD ),
					  compare_summaries_by_code    );

  return entry;

}

/*  MOD021  */
/*  MOD022  */
/********************************************************************************/
/*                  calc_sdi_max                                                */
/********************************************************************************/
/*  Description :   calculate maximum sdi for each plot                         */   
/*  Author      :   Martin W. Ritchie & Jeff D. Hamann                          */
/*  Date        :   January 4, 2000                                             */
/*  Returns     :   void                                                        */
/*  Comments    :                                                               */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code      - return code for calling function to check */
/*  vector<double> *coeffs_ptr      -   pointer to a vector of doubles that     */
/*                                      contain the coefficients for the        */
/*                                      functional species code                 */
/********************************************************************************/
/*  Formula :                                                                   */
/*  Source  : SYSTUM-1                                                          */
/*  Coeffs  :                                                                   */
/********************************************************************************/

void __stdcall calc_max_sdi(
			    unsigned long   *return_code,
			    unsigned long           n_species,
			    struct SPECIES_RECORD   *species_ptr,
			    unsigned long           n_coeffs,
			    struct COEFFS_RECORD    *coeffs_ptr,
			    unsigned long           n_plants,
			    struct PLANT_RECORD     *plants_ptr,
			    unsigned long           n_points,
			    double                  *max_sdi  )
{

  /* local variables */
  unsigned long           i;
  double                  wtd_sum_sdimx;
  double                  sum_ba;

  unsigned long           n_sp_sums;
  struct COEFFS_RECORD    *c_ptr = NULL;
  struct SUMMARY_RECORD   *temp_sp_ptr = NULL;
  struct SUMMARY_RECORD   *sp_sum_ptr = NULL;

  wtd_sum_sdimx   = 0.0;
  sum_ba          = 0.0;
  *max_sdi        = 0.0;
  *return_code = CONIFERS_SUCCESS; /* this wasn't getting set */

  /* build an array of the functional species summaries */
  /* testing section for calculating the fsp_summaries */
  //Rprintf( "before build_species_summaries(rc=%d) @ %s, %d\n", *return_code, __FILE__, __LINE__ );
  temp_sp_ptr = build_species_summaries(   return_code, 
					   n_species, 
					   species_ptr,
					   n_plants,
					   plants_ptr,
					   &n_sp_sums );

  //Rprintf( "after build_species_summaries(rc=%d) @ %s, %d\n", *return_code, __FILE__, __LINE__ );

  /* need_error_trap_here */
  /* check the return code here */
  if( *return_code != CONIFERS_SUCCESS )
    {
      /* print some heinous message */
      //Rprintf( "build_species_summaries return_code != CONIFERS_SUCCESS %s, %d\n", __FILE__, __LINE__ );
    }



  update_species_summaries(return_code, 
			   n_species,
			   species_ptr,
			   n_coeffs,
			   coeffs_ptr,
			   n_plants,
			   plants_ptr,
			   n_points,
			   n_sp_sums,
			   temp_sp_ptr );
  /* need_error_trap_here */
  /* check the return code here */
  if( *return_code != CONIFERS_SUCCESS )
    {
      /* print some heinous message */
      //Rprintf( "update_species_summaries return_code != CONIFERS_SUCCESS %s, %d\n", __FILE__, __LINE__ );
    }

  sp_sum_ptr = &temp_sp_ptr[0];
  for( i = 0; i < n_sp_sums; i++, sp_sum_ptr++ )
    {
      /* get the coeff entry for the current coeffs   */
      /* and check the type for c,h,s,etc and also    */
      /* get the max sdi value for the entry          */
      c_ptr = &coeffs_ptr[species_ptr[sp_sum_ptr->code].fsp_idx];

      /* MOD031 */
      /* don't include non-stocked functional species */
      if( is_tree( c_ptr ) )
	{
	  wtd_sum_sdimx +=    sp_sum_ptr->basal_area * ( log( 10.0 ) + 
							 log( species_ptr[sp_sum_ptr->code].max_sdi ) / REINEKE_B1 );
	  sum_ba        += sp_sum_ptr->basal_area;
	}
    } 

   
  /* print out some diagnostics here */
   

  /* print some heinous message */
  //Rprintf( "sum_ba = %lf\n", sum_ba );
   
  if( sum_ba > 0.0 )
    {
      *max_sdi = exp( ( ( wtd_sum_sdimx / sum_ba ) - log( 10.0 ) ) * REINEKE_B1 );
    }
  else
    {
      *max_sdi = 450.0;
    }

  free( temp_sp_ptr );

}


/********************************************************************************/
/*                  calc_sites_from_awi                                         */
/********************************************************************************/
/*  Description :   calculates site indecies from awi values                    */   
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   September 1, 2000                                           */
/*  Returns     :   void                                                        */
/*  Comments    :                                                               */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  double awi - water holding capacity * log( precip )                         */
/*  double *df_site - pointer to variable that stores douglas-fir site index    */
/*  double *pp_site - pointer to variable that stores p. pine    site index     */
/*  THIS FUNCTION IS NO LONGER USED                                             */
/*  IF THERE ARE ANY CALLS TO THIS FUNCTION, WE HAVE PROBLEMS                   */
/********************************************************************************/
/*  Formula :                                                                   */
/*  Source  : NA                                                                */
/*  Coeffs  :                                                                   */
/********************************************************************************/
void calc_sites_from_awi( 
			 double  awi, 
			 double  *df_site, 
			 double  *pp_site )
{

  *df_site = 0.0;
  *pp_site = 0.0;

  *df_site =  259.7758177637928 - 
    14.4264837081033 * 
    sqrt( 165.5600588801375 - 1.0 * awi );

  *pp_site =  199.0372872872872 - 
    12.9164043048688 * 
    sqrt( 96.1143063337087 - 1.0 * awi );

}

/*  MOD038  */

/********************************************************************************/
/*                  calc_sites_from_whc                                         */
/********************************************************************************/
/*  Description :   calculates site indecies from awi values                    */   
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   September 1, 2000                                           */
/*  Returns     :   void                                                        */
/*  Comments    :                                                               */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long - *return_code                                                */
/*  double whc - water holding capacity; in inches                              */
/*  double precip - annual precip in inches                                     */
/*  double *df_site - pointer to variable that stores douglas-fir site index    */
/*  double *pp_site - pointer to variable that stores p. pine    site index     */
/********************************************************************************/
/*  Formula : segmented linear                                                  */
/*  Source  : SI_temp_fits.xls                                                  */
/*  Coeffs  : imbedded                                                          */
/********************************************************************************/
void calc_sites_from_whc(
			 unsigned long *return_code, 
			 double  whc,
			 double  precip, 
			 double  *df_site, 
			 double  *pp_site )
{
  *df_site = 0.0;
  *pp_site = 0.0;

  /* MOD040 - rewrote the whole function, mwr*/

  *return_code = CONIFERS_SUCCESS;
  /* this shouldnt happen because values > 15 are not accepted  */
  if(whc > 16.0)
    {
      *df_site = 150.0;
      *pp_site = 150.0;
      return;
    }
  /*  if whc is less than 1.8, this is unacceptable value of whc  */
  if ( whc < 1.8 )
    {
      *return_code = FILL_SAMPLE_FAILED;
      return;    
    }

  if ( whc > 5.0 )
    {
      *df_site =  56.9 + 6.00 * whc;
      *pp_site =  56.9 + 6.00 * whc;
    }
  else          /* for values of whc between 1.8 and 5 */
    {
      *df_site =  27.11 + 12.05 * whc;
      *pp_site =  27.0 + 12.05 * whc;
    } 

  /* for now pp_site is the same as df_site, this should change */
  /* this temporary fix is in SI_temp_fits.xls                  */
  return;

}



/********************************************************************************/
/*                  get_type_count                                              */
/********************************************************************************/
/*  Description :   function returns the number of observations in the plant    */   
/*                  list for the given type                                     */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   September 1, 2000                                           */
/*  Returns     :   void                                                        */
/*  Comments    :                                                               */
/*  Arguments   :                                                               */
/*  return count    - number of observations in the plant list that match a type*/
/********************************************************************************/
/*  Formula :                                                                   */
/*  Source  : NA                                                                */
/*  Coeffs  :                                                                   */
/********************************************************************************/
void get_type_count( 
		    unsigned long           *return_code,
		    unsigned long           type,
		    unsigned long           *type_count,
		    unsigned long           n_species,
		    struct SPECIES_RECORD   *species_ptr,
		    unsigned long           n_coeffs,
		    struct COEFFS_RECORD    *coeffs_ptr,
		    unsigned long           n_plants,
		    struct PLANT_RECORD     *plants_ptr )
{

  unsigned long           i;
  struct  PLANT_RECORD    *plant_ptr;
  struct  COEFFS_RECORD   *c_ptr;

  *type_count = 0L;

  /* start filling in the array from the plant list               */
  /* first find the number of unique species in the plant list    */
  plant_ptr = &plants_ptr[0];
  for( i = 0; i < n_plants; i++, plant_ptr++ )
    {

      c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

      if( c_ptr->type == type )
	{
	  (*type_count)++;
	}

    }

}



/********************************************************************************/
/*   get_in_taller_attribs                                                      */
/********************************************************************************/
/*  Description :   function returns the vectors, by plant type for the two     */   
/*                  attributes (bait/cait) in taller                            */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   December 19, 2007                                           */
/*  Returns     :   void                                                        */
/*  Comments    :   the bins are scaled linearly between zero and  MAX_AIT_THT  */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  struct PLANT_RECORD     *plants_ptr a pointer to a plant structure          */
/*  struct PLOT_RECORD      *plot_ptr,                                          */
/*  double  bait[PLANT_TYPES][CON_BAT_SIZE];  basal area in taller              */
/*  double  cait[PLANT_TYPES][CON_BAT_SIZE];  crown area in taller              */
/********************************************************************************/
/*  Formula :                                                                   */
/*  Source  : NA                                                                */
/*  Coeffs  :                                                                   */
/********************************************************************************/
void get_in_taller_attribs(
			   struct PLANT_RECORD     *plant_ptr,
			   struct PLOT_RECORD      *plot_ptr,
			   double                  *bait,
			   double                  *cait )
{

  unsigned long i;
  unsigned long htidx;

  /* null out the output vectors */
  memset( bait, 0, sizeof( double ) * PLANT_TYPES );
  memset( cait, 0, sizeof( double ) * PLANT_TYPES );

  /* get the slope to scale the tht into the array */
  htidx = (int)( (plant_ptr->tht-0.5) / ((double)AIT_BIN_RES));

  if( htidx < 0 )
    {
      htidx = 0;
    }

  if( htidx >= AIT_SIZE )
    {
      htidx = AIT_SIZE - 1;
    }

  for( i = 0; i < PLANT_TYPES; i++ )
    {
      bait[i] = plot_ptr->bait[i][htidx];
      cait[i] = plot_ptr->cait[i][htidx];
    }

}

/********************************************************************************/
/*   get_in_larger_attribs                                                      */
/********************************************************************************/
/*  Description :   function returns the vectors, by plant type for bal         */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   August 1, 2014                                              */
/*  Returns     :   void                                                        */
/*  Comments    :   the bins are scaled linearly between zero and  MAX_AIT_THT  */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  struct PLANT_RECORD     *plants_ptr a pointer to a plant structure          */
/*  struct PLOT_RECORD      *plot_ptr,                                          */
/*  double  bal[PLANT_TYPES][CON_BAT_SIZE];  basal area in larger               */
/********************************************************************************/
/*  Formula :                                                                   */
/*  Source  : NA                                                                */
/*  Coeffs  :                                                                   */
/********************************************************************************/
void get_in_larger_attribs(
			   struct PLANT_RECORD     *plant_ptr,
			   struct PLOT_RECORD      *plot_ptr,
			   double                  *bal )
{

  unsigned long i;
  unsigned long htidx;

  /* null out the output vectors */
  memset( bal, 0, sizeof( double ) * PLANT_TYPES );

  /* get the slope to scale the tht into the array */
  htidx = (int)( (plant_ptr->dbh) / ((double)AIT_BIN_RES));

  if( htidx < 0 )
    {
      htidx = 0;
    }

  if( htidx >= AIT_SIZE )
    {
      htidx = AIT_SIZE - 1;
    }

  for( i = 0; i < PLANT_TYPES; i++ )
    {
      bal[i] = plot_ptr->bal[i][htidx];
    }

}


void set_in_taller_attribs(
			   struct PLANT_RECORD     *plant_ptr,
			   struct COEFFS_RECORD    *c_ptr,
			   struct PLOT_RECORD      *plot_ptr )
{

  unsigned long k;
  unsigned long htidx;

  /* compute (fill in) the bal and cait arrays */
  htidx = (int)( (plant_ptr->tht-0.5) / ((double)AIT_BIN_RES));


  /* limit the array to a 300 feet */
  if( htidx < 0 )
    {
      htidx = 0;
    }

  if( htidx > AIT_SIZE )
    {
      htidx = AIT_SIZE ;
    }

  /* fill in the "values in taller" array */
    /* where the basal area in taller values are above 6 inches */
    for( k = 0; k < htidx; k++ )
      {
	plot_ptr->bait[c_ptr->type][k] += ( plant_ptr->d6 * plant_ptr->d6 * FC_I * plant_ptr->expf );
	plot_ptr->cait[c_ptr->type][k] += ( plant_ptr->crown_area * plant_ptr->expf );
      }

}
/* this is the code for basal area in larger for CIPS*/
/* made by MWR August 1, 2014                */
void set_in_larger_attribs(
                           struct PLANT_RECORD    *plant_ptr,
                           struct COEFFS_RECORD   *c_ptr,
                           struct PLOT_RECORD     *plot_ptr)
{

  unsigned long k;
  unsigned long htidx;

  /* compute (fill in) the bal  array */
  htidx = (int)( (plant_ptr->dbh) / ((double)AIT_BIN_RES));


  /* limit the array to a 75 inches feet */
  if( htidx < 0 )
    {
      htidx = 0;
    }

  if( htidx > AIT_SIZE )
    {
      htidx = AIT_SIZE ;
    }

  /* fill in the "values in taller" array */
    /* where the basal area in taller values are above 6 inches */
    for( k = 0; k < htidx; k++ )
      {
	    plot_ptr->bal[c_ptr->type][k] += ( plant_ptr->dbh * plant_ptr->dbh * FC_I * plant_ptr->expf );
      }


}


/********************************************************************************/
/* impute                                                                       */
/********************************************************************************/
/*  Description :   imputes plant values depending on the variant's imputation  */
/*                  function, which is defined by the variant author            */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   November 21, 2011                                           */
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
void __stdcall impute_missing_values( 
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
      double                  min_prism_dbh,
	  double                  baf)
{

      /* each variant has it's own impute function  */
      switch( variant )
      {
	    case CONIFERS_SWO:
            /* this function is in model_swo.c */
            swo_impute( 
			        return_code,
			        n_species,
			        species_ptr,
			        n_coeffs,
			        coeffs_ptr,
			        variant,
			        n_plants,
			        plants_ptr,
			        n_points,
			        plots_ptr,
			        fixed_plot_radius,
			        min_prism_dbh,
			        baf );
            break;

	    case CONIFERS_SMC:
            /* this function is in model_smc.c */
            smc_impute( 
			        return_code,
			        n_species,
			        species_ptr,
			        n_coeffs,
			        coeffs_ptr,
			        variant,
			        n_plants,
			        plants_ptr,
			        n_points,
			        plots_ptr,
			        fixed_plot_radius,
			        min_prism_dbh,
			        baf );
            break;

	    case CONIFERS_SWOHYBRID:
            /* this function is in model_swohybrid.c */
            swohybrid_impute( 
			        return_code,
			        n_species,
			        species_ptr,
			        n_coeffs,
			        coeffs_ptr,
			        variant,
			        n_plants,
			        plants_ptr,
			        n_points,
			        plots_ptr,
			        fixed_plot_radius,
			        min_prism_dbh,
			        baf );
            break;

	    case CONIFERS_CIPS:
            /* this function is in model_cips.c */
            cips_impute( 
			        return_code,
			        n_species,
			        species_ptr,
			        n_coeffs,
			        coeffs_ptr,
			        variant,
			        n_plants,
			        plants_ptr,
			        n_points,
			        plots_ptr,
			        fixed_plot_radius,
			        min_prism_dbh,
			        baf );
            break;

	    default:
            /* this function is in model_swo.c */
            swo_impute( 
			        return_code,
			        n_species,
			        species_ptr,
			        n_coeffs,
			        coeffs_ptr,
			        variant,
			        n_plants,
			        plants_ptr,
			        n_points,
			        plots_ptr,
			        fixed_plot_radius,
			        min_prism_dbh,
			        baf );
      }

  if( *return_code != CONIFERS_SUCCESS )
  {
	*return_code = FILL_VALUES_ERROR;  
	return;
  }


}



/* this interface is deprecated and has been                                */
/* replaced/renamed inpute_missing_values()                                 */
void __stdcall fill_in_missing_values( 
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
      double                  min_prism_dbh,
	  double                  baf)
{

    impute_missing_values( 
      return_code,
      n_species,
      species_ptr,
      n_coeffs,
      coeffs_ptr,
      variant,
      n_plants,
      plants_ptr,
      n_points,
      plots_ptr,
      fixed_plot_radius,
      min_prism_dbh,
	  baf);

}

/* moved from the variant specific file(s) - moved from the swo_model.c */


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
void calc_biomass(
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
void calc_volume(
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


