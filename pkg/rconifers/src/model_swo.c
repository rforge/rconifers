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

/* 	$Id: model_swo.c 930 2014-01-29 21:48:28Z mritchie $	 */

/* #ifndef lint */
/* static char vcid[] = "$Id: model_swo.c 930 2014-01-29 21:48:28Z mritchie $"; */
/* #endif /\* lint *\/ */




//#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"


/* diameter calculations    */
/* S1 crown width equations */
static void swo_calc_crown_width( 
      unsigned long   *return_code,
      double          d6_area,
      double          total_height,
      double          *pred_crown_width,
      double          *pred_crown_area,
      double          *coeffs_ptr);

/* S2 coefficients: MW */
static   void calc_max_crown_width( 
      unsigned long   *return_code,
      double          dbh,
      double          total_height,
      double          *pred_max_crown_width,
      double          *coeffs_ptr);

/* S3 crown ratio calculations */
static   void swo_calc_crown_ratio( 
      unsigned long   *return_code, 
      double          total_height,
      double          d6,
      double          *pred_cr,
      double          *coeffs_ptr );


/* S4, coefficients: MS     */
static   void calc_d6_from_total_height(       
      unsigned long   *return_code,
      double          total_height, 
      double          *pred_d6,
      double          *coeffs_ptr );

static   void calc_d6_from_ht_and_dbh(       
      unsigned long   *return_code,
      double          total_height,
      double          dbh,
      double          *pred_d6,
      double          *coeffs_ptr );

//static   void calc_dbh_from_height(       
//      unsigned long   *return_code,
//      double          total_height,
//      double          *pred_dbh,
//      double          *coeffs_ptr );

 static  void calc_dbh_from_height_and_d6(       
      unsigned long   *return_code,
      double          d6,
      double          total_height,
      double          *pred_dbh,
      double          *coeffs_ptr );

 static  void calc_exp_from_cover_and_ca(
      unsigned long   *return_code,
      double          pct_cover,
      double          crown_area,
      double          *pred_expf);

/* these functions are not referenced in this model??? */

/* plant height calculations    */
// static  void calc_height_from_d6(       
//      unsigned long   *return_code,
//      double          d6,
//      double          *pred_total_height,
//      double          *coeffs_ptr );

// static  void calc_height_from_dbh(       
//      unsigned long   *return_code,
//      double          dbh,
//      double          *pred_total_height,
//      double          *coeffs_ptr );

//static   void calc_nstems_from_cw_and_height(
//      unsigned long   *return_code,
//      double          total_height,
//      double          crown_width,
//      long            *pred_n_stems,
//      double          *coeffs_ptr);


static   void calc_height_growth(
      unsigned long *return_code,	
      double        total_height, 
      double        crown_ratio,
      double        d6_area,
      double        precip, 
      double        h20_holding_capacity,
      double        slope,
      double        aspect,
      double        cacon,
      double        catcon,
      double        cahw,
      double        cash,
      double        cathw,
      double        catsh,
      double        basal_d,
      double        elevation,
      double        random_norm_0_1,
      double        random_unif_0_1a,
      double        random_unif_0_1b,
      long          ind_random,
      double        prob_browse,
      double        prob_top_damage,
      unsigned long use_precip,     
      double        *height_growth,
      double        *coeffs_ptr,
	  unsigned long plant_type);

static   void calc_d6_growth( 
      unsigned long *return_code,
      double        height_growth,
      double        crown_width,
      double        total_height,
      double        d6,
      double        cat_shrubs,
      double        cat_conifers,
      double        cat_hardwoods,
      double        ca_shrubs,
      double        ca_conifers,
      double        ca_hardwoods,
      double        h20_holding_capacity,
      double        *pred_d6_growth,
      double        *coeffs_ptr,
	  unsigned long plant_type);

 static  void calc_dbh_growth(   
      unsigned long   *return_code,
      double          total_height,
      double          height_growth,
      double          crown_ratio,
      double          current_dbh,
      double          current_d6,
      double          d6_growth,
      double          *pred_dbh_growth,
      double          *coeffs_ptr );

static   void calc_cr_growth(
      unsigned long   *return_code,
      int             cr_flag,
      double          total_height,
      double          height_growth,
      double          crown_ratio,
      double          conifer_ca,
      double          hardwood_ca,
      double          shrub_ca,
      double          unif_0_1,
      double          *cr_growth,
      double          *coeffs_ptr );

static   void calc_cw_growth(   
      unsigned long     *return_code,
      double            total_height,
      double			height_growth,
      double			crown_width,
      double			ca_conifers,
      double			ca_hardwoods,
      double			ca_shrubs,
      double			catcon,
      double		    unif_0_1,
      double            *pred_cw_growth,
      double            *coeffs_ptr,
      unsigned long     plant_type);

static void swo_calc_endemic_mortality(   
      unsigned long   *return_code,
      double          expansion_factor,
      double          *pred_mortality,
      double          *coeffs_ptr );



/********************************************************************************/
/* swo_impute                                                                   */
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
void swo_impute( 
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

  double  bait[PLANT_TYPES];
  double  cait[PLANT_TYPES];

  //double        cait_c;  //unused removed jan 2014 mwr;
  //double        cait_h;  //unused removed jan 2014 mwr;
  //double        cait_s;  //unused removed jan 2014 mwr;

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
			
			      calc_d6_from_ht_and_dbh(return_code,
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

			  
			      calc_d6_from_total_height(  return_code, 
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
		//if( plant_ptr->d12 == 0.0 && plant_ptr->tht > ( 30.0 * CM2FT ) )
		//{
        //            plot_ptr = get_plot( plant_ptr->plot, n_points, plots_ptr );
        //            cips_calc_d12_from_ht_and_veg_cov(  return_code,
		//	                    plant_ptr->tht,
		//	                    plot_ptr->shrub_pct_cover,
		//	                    &plant_ptr->d12,
		//	                    c_ptr->dbh_ht_veg_cov );
        //}

          /* this is the code for computing the missing dbh when the tree */
          /* is taller than 4.5 feet                                        */
		if( plant_ptr->d6 > 0.0 && plant_ptr->dbh == 0.0 && plant_ptr->tht > 4.5)
		  {
		    *return_code = CONIFERS_SUCCESS;
		    
			  calc_dbh_from_height_and_d6(return_code,
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
		      
			
			    calc_d6_from_total_height(  return_code, 
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
		
		  
		      /* original variant */
		      swo_calc_crown_width(   return_code,
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
		
		
		    calc_exp_from_cover_and_ca( return_code,
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
            /* maybe replace these select-case statements with function calls? */
		    get_in_taller_attribs( plant_ptr, plot_ptr, bait, cait );
            

		    //cait_c       =   cait[CONIFER];     // unused removed jan 2014;
		    //cait_h       =   cait[HARDWOOD];    // unused removed jan 2014;
		    //cait_s       =   cait[SHRUB];       // unused removed jan 2014;
		    *return_code = CONIFERS_SUCCESS;

		          swo_calc_crown_ratio(return_code,  
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
	      
		    calc_max_crown_width(  return_code,
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
/* swo_project_plant                                                            */
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
void swo_project_plant(  
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
   //double                  awi;        // cruft mwr jan 2014;

   /* tree level stats */
   //double  bat_total;   // cruft mwr jan 2014;
   //double  bat_c;
   //double  bat_h;
   //double  bat_s;
   //double  bat_c_h;    // cruft mwr jan 2014;
   double  cat_c;
   double  cat_h;
   double  cat_s;
   //double  new_d6_area;  // cruft mwr jan 2014;
   //double  new_d12_area;  // cruft mwr jan 2014;


/*    unsigned long htidx = 0; */

   double  normal;
   double  browse_random_unif_0_1;
   double  top_dam_random_unif_0_1;

   double  bait[PLANT_TYPES];
   double  cait[PLANT_TYPES];
    
//   Rprintf( "%s, %d\n", __FILE__, __LINE__ );

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
   //awi         =   plot_ptr->water_capacity * log( plot_ptr->mean_annual_precip );


    /* this function takes a vector [PLANT_TYPES][TALLER_PLANT_SIZE] */
    //get_taller_attribs( plant_ptr->tht, 
    //                    plot_ptr, 
    //                    bait, 
    //                    cait );
   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );

    get_in_taller_attribs( plant_ptr, 
                            plot_ptr, 
                            bait, 
                            cait );

   //bat_c       =   bait[CONIFER];  //cruft mwr jan 2014;
   //bat_h       =   bait[HARDWOOD];
   //bat_s       =   bait[SHRUB];

   cat_c       =   cait[CONIFER];
   cat_h       =   cait[HARDWOOD];
   cat_s       =   cait[SHRUB];

   // bat_c_h     =   bat_c + bat_h; // cruft jan 2014 mwr;
   // bat_total   =   bat_c + bat_h + bat_s;  // cruft jan 2014 mwr;

   plant_ptr->tht_growth = 0.0;
   plant_ptr->d6_growth  = 0.0;
   plant_ptr->d12_growth  = 0.0;

   plant_ptr->cw_growth  = 0.0;
   plant_ptr->cr_growth  = 0.0;
   plant_ptr->dbh_growth = 0.0;
   plant_ptr->expf_change= 0.0;

//Rprintf( "%s, %d\n", __FILE__, __LINE__ );
	if( is_tree( c_ptr ) || is_shrub( c_ptr))
	{
		calc_height_growth( return_code,
			                plant_ptr->tht, 
			                plant_ptr->cr,
			                plant_ptr->d6_area,
			                plot_ptr->mean_annual_precip, 
			                plot_ptr->water_capacity,
			                plot_ptr->slope,
			                plot_ptr->aspect,
						
			                plot_ptr->ca_c,         /* crown area in conifers           */
			                cat_c,                  /* crown area in taller conifers    */

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
		calc_d6_growth(	return_code,
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
      calc_dbh_growth(  return_code,
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
      calc_cr_growth(   return_code,
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

      calc_max_crown_width( return_code,
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
   // new_d6_area = plant_ptr->d6 * plant_ptr->d6 * FC_I; // cruft mwr jan 2014
   
   if( is_tree( c_ptr ) || is_shrub( c_ptr))
   {
		calc_cw_growth( return_code,
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

//Rprintf( "before endemic_mortality == 1 check %s, %d\n", __FILE__, __LINE__ );
   if( endemic_mortality == 1 )    
   {

//Rprintf( "before swo_calc_endemic_mortality %s, %d\n", __FILE__, __LINE__ );


      swo_calc_endemic_mortality(	return_code,
      plant_ptr->expf,
      &plant_ptr->expf_change,
      &species_ptr[plant_ptr->sp_idx].endemic_mortality);

//Rprintf( "after swo_calc_endemic_mortality %s, %d\n", __FILE__, __LINE__ );
      
      if( *return_code != CONIFERS_SUCCESS )
      {
//Rprintf( "error in swo_calc_endemic_mortality %s, %d\n", __FILE__, __LINE__ );
	 return;
      }
   }
//Rprintf( "%s, %d\n", __FILE__, __LINE__ );

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
void swo_calc_crown_width( 
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
            /* MOD011 */
/*            b1 = 1.0 - exp(total_height * b1);     */
/*            new_crown_area =b0 * b1 * pow( d6_area, b2 );   */

            /* MOD014  */
            *pred_crown_area = exp( b0 + 
                                    b1 * log( d6_area * 144.0 ) + 
                                    b2 * log( total_height ) );
		if(*pred_crown_area > 2827.0)
		{
			*pred_crown_area = 2827.0; /* limit crown width to 60 feet */
		}
            /* MOD003 */
            *pred_crown_width = sqrt( *pred_crown_area * ONE_OVER_PI * 4.0); 
        }

        *return_code = CONIFERS_SUCCESS;

        /* MOD006   */
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

void calc_max_crown_width( 
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
void swo_calc_crown_ratio(
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
void calc_d6_from_total_height(       
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
void calc_d6_from_ht_and_dbh(       
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
void calc_dbh_from_height(       
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
void calc_dbh_from_height_and_d6(       
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
void calc_exp_from_cover_and_ca(
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
void calc_height_from_d6(       
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
void calc_height_from_dbh(       
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
void calc_nstems_from_cw_and_height(
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
/*                  calc_d6_growth          D1                                  */
/********************************************************************************/
/*  Description :   calc_d6_growth                                              */    
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   March 17, 2004                                              */
/*  Returns     :   void                                                        */
/*  Comments    :   refit by martin                                             */
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
/*  Coeffs  : D3                                                                */
/********************************************************************************/
void calc_d6_growth(
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
    //double  b8;
    //double  b9;
    //double  b10;
    //double  b11;
    //double  b12;
    //double  b13;
    //double  b14;
    //double  b15;

    double  c0;
    double  c1;
    double  c2;
    double  c3;
    double  c4;

    double  temp_cat_total;
    double  sum_cat;
    //double  cat_hardwoods_conifers;  // cruft mwr jan 2014;
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
    //b8  = 0.0;    // cruft mwr Jan 2014;
    //b9  = 0.0;
    //b10 = 0.0;
    //b11 = 0.0;
    //b12 = 0.0;
    //b13 = 0.0;
    //b14 = 0.0;
    //b15 = 0.0;
    c0  = 0.0;
    c1  = 0.0;
    c2  = 0.0;
    c3  = 0.0;
    c4  = 0.0;


    sum_cat                = cat_shrubs + cat_conifers + cat_hardwoods; 
    //cat_hardwoods_conifers = cat_conifers + cat_hardwoods; // cruft mwr 2014;

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
	    b0  = coeffs_ptr[0];
		b1  = coeffs_ptr[1];
		b2  = coeffs_ptr[2];
		b3  = coeffs_ptr[3];
		b4  = coeffs_ptr[4];
		b5  = coeffs_ptr[5];
		b6  = coeffs_ptr[6];
		b7  = coeffs_ptr[7];
		//b8  = coeffs_ptr[8];
		//b9  = coeffs_ptr[9];
		//b10 = coeffs_ptr[10];
		//b11 = coeffs_ptr[11];
		//b12 = coeffs_ptr[12];
		//b13 = coeffs_ptr[13];
		//b14 = coeffs_ptr[14];
		//b15 = coeffs_ptr[15];

		dg_trees =   pow(height_growth, b1)
					*exp(   b0                           
							+	b2 * sqrt(crown_width)      
							+	b3 * d6                     
							+	b4 * (cat_conifers/SQ_FT_PER_ACRE)*(cat_conifers/SQ_FT_PER_ACRE)    
							+	b5 * (cat_hardwoods/SQ_FT_PER_ACRE)*(cat_hardwoods/SQ_FT_PER_ACRE)  
							+	b6 * (cat_shrubs/SQ_FT_PER_ACRE)*(cat_shrubs/SQ_FT_PER_ACRE)        
							+	b7 * d6*d6);

		*pred_d6_growth = dg_trees ;
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

/*  MOD025   */
/********************************************************************************/
/*                  calc_dbh_growth        D2                                   */
/********************************************************************************/
/*  Description :   calc_dbh_growth                                             */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double         total_height     - total tree height                         */
/*  double         height_growth    - predicted height increment (annual)       */
/*  double         crown_ratio      - crown ratio (currently not used)          */
/*  double         current_dbh      - initial dbh (inches)                      */
/*  double         current_d6       - initial d6 (inches)                       */
/*  double         d6_growth        - predicted d6 growth (inches)              */
/*  double         *pred_dbh_growth - predicted dbh growth (inches)             */
/*  vector<double> *coeffs_ptr      -   pointer to a vector of doubles that     */
/*                                      contain the coefficients for the        */
/*                                      functional species code                 */
/********************************************************************************/
/*  Formula : dbhgro = d6gro * (b0 + exp(b1 + b2 * height))                     */ 
/*  Source  :                                                                   */
/*  Coeffs  : static  coeffs number 7 d6=f( h dbh)                              */
/********************************************************************************/
void calc_dbh_growth(   
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
void calc_height_growth(
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
    double          *height_growth,
    double          *coeffs_ptr,
	unsigned long   plant_type)
{

/****************************  right now b25 and b26 are spares ******************************************/
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
    // double  b14; // cruft mwr jan 2014;
    // double  b15; // cruft mwr jan 2014;

    double  whc;
	double  height_var;
	double  height_for_error;
    int     broken;
    int     browsing;
    double  cat;
    double  temp_growth;

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

    b0                  = 0.0;
    b1                  = 0.0;
    b2                  = 0.0;
    b3                  = 0.0;
    b4                  = 0.0;
    b5                  = 0.0;
    b6                  = 0.0;
    b7                  = 0.0;
    b8                  = 0.0;
    b9                  = 0.0;
    b10                 = 0.0;
    b11                 = 0.0;
    b12                 = 0.0;
    b13                 = 0.0;
    //b14                 = 0.0; // cruft mwr jan 2014;
    //b15                 = 0.0; // cruft mwr jan 2014;
    
    broken              = 0;
    browsing            = 0;

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
    if (plant_type == CONIFER || plant_type == HARDWOOD) 
    {
		b0  = coeffs_ptr[0];
		b1  = coeffs_ptr[1]  * whc;
		b2  = coeffs_ptr[2]  * log(precip);                        /* Not in the model now b2 should=0     */
		b3  = coeffs_ptr[3]  * log (total_height);
		b4  = coeffs_ptr[4]  * pow(total_height, 1.50);                
		b5  = coeffs_ptr[5]  * log(crown_ratio);                       
		b6  = coeffs_ptr[6]  * (catcon/SQ_FT_PER_ACRE)*(catcon/SQ_FT_PER_ACRE);
		b7  = coeffs_ptr[7]  *  (cahw/SQ_FT_PER_ACRE) * (cahw/SQ_FT_PER_ACRE);
		b8  = coeffs_ptr[8]  *  (cash/SQ_FT_PER_ACRE) * (cash/SQ_FT_PER_ACRE);
		b9  = coeffs_ptr[9]  * (sqrt(height_for_error/2.0) * height_var);
		b10 = (double)coeffs_ptr[10] * (double)broken;
		b11 = (double)coeffs_ptr[11] * (double)browsing;
		b12 = coeffs_ptr[12] * (1.0 / total_height)  ;                         
		b13 = coeffs_ptr[13] * height_var / sqrt(2.0);             /*  alternate constant error structure   */
                                                                   /* currently unused here                 */

        temp_growth =  b12 + exp(b0 + b1 + b3 + b4 + b5 + b6 + b7 + b8) + (b9 + b13);

        if( (temp_growth) <= 0.0 )
        {
            *height_growth  = 0.0 + (b10 + b11);
        }
        else
        {
	        *height_growth  = temp_growth + (b10 + b11);          /* damage adjustments                    */
        }
    }
/******************** this concludes the tree height growth model  ******************************************/
    else if(plant_type == SHRUB )        
    {
		b0  = coeffs_ptr[0];
		b1  = coeffs_ptr[1];
		b2  = coeffs_ptr[2];
		b3  = coeffs_ptr[3];
		b4  = coeffs_ptr[4];

		*height_growth= b0/total_height 
                  + exp(b1 + b2*(log(total_height)) + b3*total_height*basal_d + b4*cat*cat);

        /****for unknown brush species set growth to zero *****/
        if (b1 == 0 && b2 == 0 && b3 == 0 && b4 == 0)
        {
            *height_growth=0.0;
        }
    }
	else
	{
		*height_growth=0.0;
	}
/******************* this concludes the shrub height growth model *******************************************/
/************************************************************************************************************/
	if( total_height + *height_growth <= 0.5)                       
        /* if predicted height is <0.5 then dont grow       */
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
void calc_cr_growth(
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
///  This needs to be checked mwr Jan 2014 pending check;
///  I checked this with the original in Jan and it looks ok;
///  I think we need to refit. This is kind of clunky.;
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
void calc_cw_growth(   
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
	double  c3;
	double  c4;
	double  c5;
	double  c6;
	//double  c7;  // cruft mwr jan 2014;
	double  c8;
	double  c9;
	double  c10;

/* this function currently does not have a random damage component but that is */
/*  what the extra parameters (c7, c8) are for and the uniform_0_1 variable    */

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
	ca_conifers = (ca_conifers/SQ_FT_PER_ACRE)*(ca_conifers/SQ_FT_PER_ACRE);
	ca_hardwoods= (ca_hardwoods/SQ_FT_PER_ACRE)*(ca_hardwoods/SQ_FT_PER_ACRE);
	ca_shrubs   = (ca_shrubs/SQ_FT_PER_ACRE)*(ca_shrubs/SQ_FT_PER_ACRE);
	catcon      = catcon/SQ_FT_PER_ACRE;
    
	c0      = coeffs_ptr[0];
	c1      = coeffs_ptr[1];
	c2      = coeffs_ptr[2];
	c3	    = coeffs_ptr[3];
	c4	    = coeffs_ptr[4];
	c5	    = coeffs_ptr[5];
	c6	    = coeffs_ptr[6];
	//c7	    = coeffs_ptr[7];  // cruft mwr jan 2014;
	c8	    = coeffs_ptr[8];
	c9      = coeffs_ptr[9];
	c10     = coeffs_ptr[10];

/*  for trees,  hw & conifers c8 is an indicator parameter 0 for tree */
	if(plant_type == CONIFER || plant_type == HARDWOOD)
	{
		temp_cwg = pow(height_growth,c1)
			      * ( c0 
				+ c2*sqrt(crown_width) 
				+ c3*ca_conifers 
				+ c4*ca_hardwoods 
				+ c5*ca_shrubs 
				+ c6*log(crown_width));

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

	else if(plant_type == SHRUB)   /* if it is a shrub, change to a different function using c8 c9 c10*/
	{
	    temp_cwg = (height_growth*2.0) * c8*pow(total_height, c9) * exp(c10*catcon);
        if(temp_cwg < 0.0)
        {
            *pred_cw_growth =0.0;
            *return_code    = CONIFERS_SUCCESS;
            return;
        }
	    *pred_cw_growth  = 0.5* temp_cwg;
	    *return_code     = CONIFERS_SUCCESS;
	    return;
	}
    else
    {
        *pred_cw_growth  = 0.0;
	    *return_code     = CONIFERS_SUCCESS;
	    return;
    }

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
void swo_calc_endemic_mortality(   
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
