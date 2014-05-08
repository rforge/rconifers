/****************************************************************************/
/*                                                                          */
/* cips_model.c                                                             */
/* functions used to predict values for the CONIFERS growth model           */
/* for the model defined in:                                                */
/* Center for Intensive Planted-Forest Silviculture. Annual Report 2010.    */
/* Oregon State University                                                  */
/* Douglas Maguire and Doug Mainwaring.                                     */
/*                                                                          */
/****************************************************************************/

/*  $Id: model_cips.c 931 2014-04-30 18:13:11Z mritchie $  */



#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"


/****************************************************************************/
/* functions in ported from jim flewelling's fortran code                   */
/****************************************************************************/
#define FLEWELLING_ET_AL_2002_STOCK_AGE 2.0f    /* stock age at planting    */
#define FLEWELLING_ET_AL_2002_STOCK_HT0 1.4f    /* stock height at planting */

/* these are macros that are only used here */
#define MAX(X,Y) ( (X) > (Y) ? (X) : (Y) )
#define BOUND(X, L, U) ( (X) < (L) ? (L) : ( (X) > (U) ? (U) : (X) ) )
#define INCDOWN(X, place) X -= pow(10.0f, -(place) )

void get_delta_h_pot (
    double    psi,          /* PSI  = PSI computed elsewhere                */
    double    ht0,          /* HT0  = avg ht of the planted stock (ft)      */
    double    htx,          /* HTx = current top height                     */
    double   *delta_output ); /* output estimate of potential 1-year growth */

static void get_gea (
    double    psi, /* PSI  = PSI computed elsewhere                         */
    double    ht0, /* HT0  = average height of the planted stock (ft)       */
    double    htx, /* HTx = current top height                              */
    double   *gea_output ); /* output growth effective age                   */

void get_psi (
   double  x,           /* x    = Number of years since planting. Make X=30-stock age   */
   double ht0,          /* HT0  = average height of the planted stock (ft)      */
   double htx,          /* HTX  = top height at age x                           */
   double *psi_output );  /* output estimate of psi                               */

static void flewelling_site_index(
		double  psi,
		double  ht0,
		double  x,
        double  *htop_output );

/****************************************************************************/
/* functions defining the model bahavior                                    */
/****************************************************************************/
static void cips_calc_height_growth(
    unsigned long           *return_code,   
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *c_ptr,
    double                  h40,    /* plot level variable */
    unsigned long           *n_years_after_planting );

static void cips_calc_dbh_growth(
    unsigned long           *return_code,   
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  ht40 );


static void cips_calc_dbh_from_ht_and_veg_cov( 
    unsigned long           *return_code,   
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr );


/* these functions:                 */
/*
static void cips_calc_d6_from_ht_and_veg_cov( 
    unsigned long           *return_code,   
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr );

static void cips_calc_d12_from_ht_and_veg_cov(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr );
*/

/* with this new function */
static void cips_calc_dob_hi( 
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  hi,
    double                  *di );

static void cips_calc_cw_growth(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr,
    unsigned long           plantation_age );

// todo: these will need to be moved into the plot record
//    double                  conifer_d12_basal_area, -- these are at the plot level */
//    double                  catcon );     /* another plot level variable */

static void cips_calc_crown_ratio(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr );

static void cips_calc_max_crown_width(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr );

static void cips_calc_crown_width(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr );


static void cips_calc_endemic_mortality(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  h40 );


static int compare_plants_by_plot_plant(
    const void *ptr1,
    const void *ptr2 );



/****************************************************************************/
/* local function definitions                                               */
/****************************************************************************/



/* the variables - or some of them
where   D   =   Dbh (mm)
    D15 =   Diameter at 15 cm above ground (mm)
    D30 =   Diameter at 30 cm above ground (mm)
    H   =   Height (cm)
    h   =   Arbitrary height on stem with 0=h=H
    Hrel    =   H/H100
    H100    =   Mean height of 100 largest trees per hectare (cm)
    Vegcov  =   Competing vegetation cover (%)
    PAGE    =   Plantation age (years)
    VGAP    =   100 - initial Vegcov
    BADF    =   Initial Douglas-fir basal area at 30 cm from ground (m2/ha)

*/

/* these are the functions you have already written */
//cips_calc_height_growth - eq. 1, page 32
//cips_calc_dbh_growth - eq 4, page 32.
//cips_calc_dbh_from_ht_and_veg_cov - eq 3, page 42
//cips_calc_d6_from_ht_and_veg_cov - where did this come from?
//cips_calc_d12_from_ht_and_veg_cov - where did this come from?
//cips_calc_endemic_mortality - eq. 2, page 38


/****************************************************************************/
/* cips_project_plant                                                       */
/****************************************************************************/
/*  Description :   computes plant growth components for CIPS variant		*/
/*  Author      :   Jeff D. Hamann                                          */
/*  Date        :   July 01, 2011											*/
/*  Returns     :   void                                                    */
/*  Comments    :   updated values are height_growth, d6_growth, dbh_growth,*/
/*                  crown_ratio_growth, crown_width_growth, and             */
/*                  max_crown_width                                         */
/*  Arguments   :   unsigned long *return_code  - pointer to a return code  */
/*                  unsigned long n_plants - total number fo plants in the  */
/*                      plants pointer array                                */
/*                  struct PLANT_RECORD *plants_ptr - array of plants in the*/
/*                      sample to be projected                              */
/*                  struct PLOT_RECORD  *plot_ptr   - pointer to the current*/
/*                      plot that is to be grown                            */
/*                  unsigned long n_species - size of the species_ptr       */
/*                  struct SPECIES_RECORD   *species_ptr - array of         */
/*                      SPECIES_RECORD's that hold species specific info    */
/*                  unsigned long   n_coeffs - sizes of the coeffs_ptr array*/
/*                  struct COEFFS_RECORD *coeffs_ptr - array of coefficients*/
/*                      that are used to project the individual plants on   */
/*                      plot.                                               */
/****************************************************************************/

/* this function needs to update the following values	*/
/* update the current tree values...					*/
//plant_ptr->d6           += plant_ptr->d6_growth;
//plant_ptr->dbh          += plant_ptr->dbh_growth;
//plant_ptr->tht          += plant_ptr->tht_growth;
//plant_ptr->cr           += plant_ptr->cr_growth;
//plant_ptr->crown_width  += plant_ptr->cw_growth;
//plant_ptr->expf         -= plant_ptr->expf_change;


/* these are the variables that are updated using the */


/* you also should turn this into pseudocode using  */
/* \usepackage{algorithm}                           */
/* if possible.                                     */


/*
\begin{algorithm}[tp!]
  \SetLine 
   \KwData{$\mathbf{x},\hat{\upxi},\mathbf{f}(\mathbf{x},\hat{\upxi}),\mathbf{g}(\mathbf{x},\hat{\upxi}),\mathbf{h}(\mathbf{x},\hat{\upxi})$}  
  \KwData{yield streams, adjacency list(matrix), constraint
    formulae, software}
  \KwResult{($\mu$+$\lambda$)-PAES Multi-Objective Evolutionary
    Algorithm \citep{1108875}}  
  $g \leftarrow 1$\;
  $\mu_{g} = \emptyset$\;
  $\lambda_{g} \leftarrow GenerateRandom ()$\;

    \If{ is\_tree( C ) \Or is_shrub( C ) } {
        CIPS\_CALC\_HEIGHT\_GROWTH(Plant,Plot,C,etc.)
    }


    \If{ is\_tree( C ) } {

        \If{ is\_tree( C ) \And P^{height} > 15.0 cm } {
            cips_calc_d6_from_ht_and_veg_cov(P^{height},P^{shrub\_cover})
        }

        \If{ is\_tree( C ) \And P^{height} > 30.0 cm } {
            cips_calc_d12_from_ht_and_veg_cov(P^{height},P^{shrub\_cover})
        }

        \If{ is\_tree( C ) \And P^{height} > 4.5 ft \And P^{DBH} == NULL } {
			cips_calc_dbh_from_ht_and_veg_cov(Plant^{height},P^{shrub\_cover})
        }

  }

  $Evaluate ( \lambda_{g})$\;
  $\mu_{g} \leftarrow UpdateParetoArchive(\lambda_{g})$\;
  \For{$g \leftarrow 2$ \KwTo $G$ } {    
    $\lambda_{g} \leftarrow SelectCandidatesFromArchive(\mu_{g-1})$\;
    $\lambda_{g} \leftarrow Mutate(\lambda_{g})$\;
    $Evaluate (\lambda_{g})$\;
    $\mu_{g} \leftarrow UpdateParetoArchive(\lambda_{g})$\;
  }    
  $ExportArchive(\mu_{G})$\;
  \caption{Pareto archiving multi-objective evolutionary algorithm.}
  \label{alg:ch4_moea_paes_algorithm}
\end{algorithm}
*/


/* calculates the height growth */
    
void cips_project_plant(  
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
   unsigned long           use_rand_err,
   struct SUMMARY_RECORD   *before_sums,
   unsigned long           use_genetic_gains,
   unsigned long           genetics_age_cut,
   unsigned long			plantation_age,
   unsigned long            *n_years_projected )
{

   struct COEFFS_RECORD    *c_ptr;

   /* tree level stats */
   // double  bat_total; //unused removed jan 2014 by mwr;
   //double  bat_c;   // unused removed jan 2014 by mwr;
   //double  bat_h;   // unused removed jan 2014 by mwr;
   //double  bat_s;   // unused removed jan 2014 by mwr;
   //double  bat_c_h;  // unused removed jan 2014 by mwr;
   ///double  cat_c;  //unused, rem jan 2014 mwr;
   //double  new_d6_area; //unused, rem jan 2014 mwr;
   //double  new_d12_area; //unused, rem jan 2014 mwr;

    /* variables added for d6_growth and d12_growth */
   double   old_d6;
   double   new_d6;

   double   old_d12;
   double   new_d12;

   double   old_cr;
   double   new_cr;

   //double  normal;  //unused, rem jan 2014 mwr;
   //double  browse_random_unif_0_1;  //unused, rem jan 2014 mwr;
   //double  top_dam_random_unif_0_1; //unused, rem jan 2014 mwr;
   //double  basal_area; //unused, rem jan 2014 mwr;
   double  h40;
   //double  tpa_con_stand;  //unused, rem jan 2014 mwr;
   //double  current_height; //unused, rem jan 2014 mwr;
   //double  new_height;  //unused, rem jan 2014 mwr;
   //double  pred_crown_ratio = 0.0; //unused, rem jan 2014 mwr;


   double  bait[PLANT_TYPES];
   double  cait[PLANT_TYPES];

    
#ifdef _DEBUG

   if( *n_years_projected == 4 && plant_ptr->plant == 35 )
   {
       //fprintf( stdout, "stop\n" );
   }

#endif


   /* get the supporting structures for the plant */
   c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

   /* don't go any further if you don't    */
   /* have a functional species code       */
   if( c_ptr == NULL )
   {
      *return_code = CONIFERS_ERROR;
      return;
   }

   /* get a uniform deviate for the browse */
   /* and one for the top damage           */
   //normal                  = (double)gauss_dev();  // unused rem mwr;
   //browse_random_unif_0_1  = uniform_0_1();    // unused rem mwr;
   //top_dam_random_unif_0_1 = uniform_0_1();    // unused rem mwr;


   /* get the variables to project the plant    */
   /* these are stand level variables           */
   //tpa_con_stand  = before_sums->con_tpa;  // unused rem mwr;
   h40		      = before_sums->height_40;
   //basal_area     = before_sums->basal_area;  // unused rem mwr;

    get_in_taller_attribs(  plant_ptr,
                            plot_ptr,
                            bait,
                            cait );

   /* basal area in [Conifers/Hardwoods/Shrubs] */
   //bat_c       =   bait[CONIFER];   // unused removed mwr jan 2014;
   //bat_h       =   bait[HARDWOOD];  // unused removed mwr jan 2014;
   //bat_s       =   bait[SHRUB];     // unused removed mwr jan 2014;

   /* crown area in [Conifers/Hardwoods/Shrubs] */
   //cat_c       =   cait[CONIFER];   // unused rem mwr;

   /* and create some temporary variables you might need for this model */
   //bat_c_h     =   bat_c + bat_h; /* ba in taller con and hw (trees) rem jan 2014 mwr */
   //bat_total   =   bat_c + bat_h + bat_s;   /* ba in tallr con_hw_sh rem jan 2014     */


	/* the CIPS variant requires the basal area of all plants at D6		*/
	/* which in this case is only douglas-fir, which is the only		*/
	/* conifer (and tree) in the model									*/

   /* these are the variables for each plant record that need to be filled in */
   /* each cycle */
    /* todo: write a function null_plant_growth_rates( &plant_ptr ); */
   plant_ptr->tht_growth = 0.0;
   plant_ptr->d6_growth  = 0.0;
   plant_ptr->d12_growth  = 0.0;    /* added for CONIFERS_CIPS */
   plant_ptr->cw_growth  = 0.0;
   plant_ptr->cr_growth  = 0.0;
   plant_ptr->dbh_growth = 0.0;
   plant_ptr->expf_change= 0.0;

   if( is_non_stocked(c_ptr) || plant_ptr->expf <= 0.0)
   {
       *return_code = CONIFERS_SUCCESS;
       return;
   }

   //current_height = plant_ptr->tht;  // unused rem mwr jan 2014;
   //new_height = 0.0;                 // unused rem mwr jan 2014;

   /* if the plant record is a tree or shrub		*/
   /* then compute the height growth for the plant	*/
   if( is_tree( c_ptr ) || is_shrub( c_ptr) )
   {
		/* compute the height growth for the individual plant   */
		/* the height growth must be returned in imperial units */
	  /* coded and validated, august 4, 2011, jdh             */
		cips_calc_height_growth(    return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    h40,
                                    n_years_projected );

		if ( *return_code != CONIFERS_SUCCESS)
        {
            return;
        }
		     
   }

    /* do these functions require the ending total  */
    /* height or the beginning total height?        */
    //if( plant_ptr->plant == 0 )
    //{
    //    //fprintf( stdout, "stop\n" );
    //}

   /* if the functional species is a tree, compute the diameter growths */
   if( is_tree( c_ptr ) )
   {

        /* for trees that will grow to over breast height */
        if( ( plant_ptr->tht + plant_ptr->tht_growth ) > ( 137.0f * CM2FT ) )
        {
            /* the plant will grow above breast height                          */
            /* but won't be above breast height at the beginning of the cycle   */
            /* so you have to predict it.                                       */
            if( plant_ptr->dbh <= 0.0 )
            {
                plant_ptr->tht += plant_ptr->tht_growth;
                cips_calc_dbh_from_ht_and_veg_cov(  return_code,
                                                plot_ptr,
                                                plant_ptr,
                                                c_ptr );
                plant_ptr->dbh_growth = plant_ptr->dbh;
                //plant_ptr->dbh = 0.0;

                plant_ptr->d12_growth = 0.0;
                old_d12 = plant_ptr->d12;
                new_d12 = 0.0;
                cips_calc_dob_hi(  return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    30.0 * CM2FT,
                                    &new_d12 );
                plant_ptr->d12_growth = new_d12 - old_d12;
                plant_ptr->d12 = old_d12;
                if( plant_ptr->d12_growth < 0.0f )
                {
                    plant_ptr->d12_growth = 0.0f;
                }

                plant_ptr->d6_growth = 0.0;
                old_d6 = plant_ptr->d6;
                new_d6 = 0.0;
                cips_calc_dob_hi(  return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    15.0 * CM2FT,
                                    &new_d6 );
                plant_ptr->d6_growth = new_d6 - old_d6;
                plant_ptr->d6 = old_d6;
                if( plant_ptr->d6_growth < 0.0f )
                {
                    plant_ptr->d6_growth = 0.0f;
                }
                
                plant_ptr->tht -= plant_ptr->tht_growth;
                plant_ptr->dbh -= plant_ptr->dbh_growth;

            }
            else
            {
                /* compute dbh, d12, and d6 growth */
                cips_calc_dbh_growth(   return_code,
                                        plot_ptr,
                                        plant_ptr,
                                        c_ptr,
                                        h40 );

                plant_ptr->tht += plant_ptr->tht_growth;
                plant_ptr->dbh += plant_ptr->dbh_growth;                
                plant_ptr->d12_growth = 0.0;
                old_d12 = plant_ptr->d12;
                new_d12 = 0.0;
                cips_calc_dob_hi(  return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    30.0 * CM2FT,
                                    &new_d12 );
                plant_ptr->d12_growth = new_d12 - old_d12;
                plant_ptr->d12 = old_d12;
                if( plant_ptr->d12_growth < 0.0f )
                {
                    plant_ptr->d12_growth = 0.0f;
                }

                plant_ptr->d6_growth = 0.0;
                old_d6 = plant_ptr->d6;
                new_d6 = 0.0;
                cips_calc_dob_hi(  return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    15.0 * CM2FT,
                                    &new_d6 );
                plant_ptr->d6_growth = new_d6 - old_d6;
                plant_ptr->d6 = old_d6;
                if( plant_ptr->d6_growth < 0.0f )
                {
                    plant_ptr->d6_growth = 0.0f;
                }
                
                plant_ptr->tht -= plant_ptr->tht_growth;
                plant_ptr->dbh -= plant_ptr->dbh_growth;                

            }

        } /* end of the if tree will grow to over breast height */
        else
        {
            plant_ptr->tht += plant_ptr->tht_growth;

            plant_ptr->d12_growth = 0.0;
            old_d12 = plant_ptr->d12;
            new_d12 = 0.0;
            cips_calc_dob_hi(  return_code,
                                plot_ptr,
                                plant_ptr,
                                c_ptr,
                                30.0 * CM2FT,
                                &new_d12 );
            plant_ptr->d12_growth = new_d12 - old_d12;
            plant_ptr->d12 = old_d12;
            if( plant_ptr->d12_growth < 0.0f )
            {
                plant_ptr->d12_growth = 0.0f;
            }

            plant_ptr->d6_growth = 0.0;
            old_d6 = plant_ptr->d6;
            new_d6 = 0.0;
            cips_calc_dob_hi(  return_code,
                                plot_ptr,
                                plant_ptr,
                                c_ptr,
                                15.0 * CM2FT,
                                &new_d6 );
            plant_ptr->d6_growth = new_d6 - old_d6;
            plant_ptr->d6 = old_d6;
            if( plant_ptr->d6_growth < 0.0f )
            {
                plant_ptr->d6_growth = 0.0f;
            }
            
            plant_ptr->tht -= plant_ptr->tht_growth;
            
        }
 
        /* calculate the new crown ratio and    */
        /* calculate the difference             */
        /* todo: since this model doesn't have a crown ratio growth function, */
        /* then assume the growth is always zero?								*/
		plant_ptr->cr_growth = 0.0;
        //pred_crown_ratio = 0.0;  // unused removed by mwr jan 2014;
        old_cr = plant_ptr->cr;
        new_cr = 0.0f;
        
        /*
        cips_calc_crown_ratio( return_code,
			                plant_ptr->tht + plant_ptr->tht_growth,
                            &pred_crown_ratio,
			                c_ptr->crown_ratio );        
        */

        /* this function is ending crown ratio  */
        /* so the function needs to include     */
        /* the height growth                    */
        plant_ptr->tht += plant_ptr->tht_growth;
        cips_calc_crown_ratio(  return_code,
                                plot_ptr,
                                plant_ptr,
                                c_ptr );
        plant_ptr->tht -= plant_ptr->tht_growth;
        new_cr = plant_ptr->cr;
        plant_ptr->cr_growth = new_cr - old_cr;

        /* restrict the crowns to always receed */
        /* todo: check this assumption/behavior */
        if( plant_ptr->cr_growth > 0.0f )
        {
            plant_ptr->cr_growth = 0.0f;
        }

        if( *return_code != CONIFERS_SUCCESS )
        {
    	   return;
        }

   } /* end of the tree species function */


   /* compute the new basal area at the base of the stem			*/
   /* you might have to track the d6 and d12 growth here instead	*/
   /* are these getting used for anything? */
   //new_d6_area	= plant_ptr->d6		* plant_ptr->d6 * FC_I;
   //new_d12_area = plant_ptr->d12	* plant_ptr->d12 * FC_I;

	/* if the functional species is either a tree or a shrub	*/
	/* then compute the crown width growth. For shrubs, the		*/
	/* is the change in percent cover for the competeting veg	*/
   if( is_tree( c_ptr ) || is_shrub( c_ptr ) || is_forb( c_ptr ) )
   {
		/* this is converted from the shrubs and trees	*/
		/* todo: is there a crown width model for the trees?	*/
		/* todo: you'll need to find out if this is the doug-fir basal area */
        /* make sure all of the arguments are getting through */

       //if( *n_years_projected == 4 && plant_ptr->plant == 35 )
       //if( plant_ptr->plant == 35 )
       //{
       //    //fprintf( stdout, "stop\n" );
       //}

        /* which variables need to get into this function finally? */
       cips_calc_cw_growth(     return_code,
                                plot_ptr,
                                plant_ptr,
                                c_ptr,

	                            // from doug's new equation
	                            plantation_age );

	                            //plot_ptr->d12ba_c,
                                //plot_ptr->cait );


        if ( *return_code != CONIFERS_SUCCESS)
        {
            return;
        }
   }



   /* finally, compute the mortality */
   if( endemic_mortality == 1 )    
   {
	   /* coded and checked - jdh, august 4, 2011	*/
	   /* against the osu spreadsheet				*/
       cips_calc_endemic_mortality( return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    h40 );

                                   // &species_ptr[plant_ptr->sp_idx].endemic_mortality );

        if ( *return_code != CONIFERS_SUCCESS)
        {
            return;
        }
   }


   *return_code = CONIFERS_SUCCESS;

}





/********************************************************************************/
/*                  cips_calc_height_growth                                     */
/********************************************************************************/
/*  Description :   calc_height_growth                                          */   
/*  Author      :   Jeff D. Hamann and Doug Maguire                             */
/*  Date        :   June 16, 2011                                               */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  flew_site                   -   Flewellings site index                      */
/*  h40                         -   height of 40 tallest trees per acre         */
/*                              -   which is the same as the 100 tallest stems  */
/*                                  per hectare                                 */
/* -----------------------------------------------------------------------------*/
/*  *height_growth              -   predicted height growth                     */
/*  vector<double> *coeffs_ptr  -   pointer to a vector of doubles that         */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/*  unsigned long  plant_type       - growth form of the plant record           */
/********************************************************************************/
/*  Formula : equation 4, page 32                                               */
/*  Source  : Center for Intensive Planted-Forest Silviculture,                 */
/*              2010 Annual Report                                              */
/*              http://www.fsl.orst.edu/cips/                                   */
/*  Coeffs  : ht_growth                                                         */
/********************************************************************************/

//		cips_calc_height_growth( return_code,
//                                    plot_ptr,
//                                    c_ptr,
//                                    plant_ptr,
//                                    h40,
//                                    n_years_projected );


//                                    plot_ptr,
//                                    c_ptr,
//                                    plant_ptr,
//                                    h40,
//                                    n_years_projected );

    //double          total_height, 
    //double          flew_site,
    //double          pct_veg_cov,    /* todo: verify shrub pct cover */
    //double          *height_growth,
    //double          *coeffs_ptr,
    //unsigned long   plant_type,

    /* these will need to be removed for the release */
    //double          *debug_rel_height,
    //double          *debug_rh_mod_hg,
    //double          *debug_cv_mod_hg,
    //unsigned long   *debug_after_first_season_flag,
    //double          *debug_delta_h_pot,
    //double   *debug_psi


void cips_calc_height_growth(
    unsigned long           *return_code,   
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeff_ptr,
    double                  h40,
    unsigned long           *n_years_after_planting )
{

    double          b1          =   0.0;
    double          b2          =   0.0;
    double          b3          =   0.0;
    double          b4          =   0.0;  
    double          b5          =   0.0;

    double          rel_height  =   0.0;
    double          rh_mod_hg   =   0.0;
    double          cv_mod_hg   =   0.0;
    unsigned long   after_first_season = 0; /* for now, just keep this as after the first */
    double          delta_h_pot =   0.0;
    double          psi         =   0.0;
    
    /********** initialize variables ***********/
    *return_code = CONIFERS_SUCCESS;

    //if( coeffs_ptr == NULL )
    if( coeff_ptr == NULL )
    {
        //*height_growth  =   0.0;
        plant_ptr->tht_growth = 0.0;
        *return_code    =   CONIFERS_ERROR;
        return;
    }

    /* todo: make sure this is the correct formula for the relateive height */
    //*height_growth      = 0.0;              /* set to zero for default      */
    plant_ptr->tht_growth = 0.0;
    rel_height           = ( plant_ptr->tht * FT2CM ) / ( h40 * FT2CM );  /* compute the relative height  */

    if( *n_years_after_planting > 0 )
    {
        after_first_season = 1;
    }

    /* I = 0 for the first season after planting, 1 otherwise */
    /* for now, assume I = 1 always */
    //first_season_ind = 1;

    /****************** check for some errors *****************/
    if( plant_ptr->tht < 0.0)
    {
        //*height_growth    = 0.0;
        plant_ptr->tht_growth = 0.0;
        *return_code        = INVALID_INPUT_VAL;
        return;
    }

    /* todo: you might want to confirm this with Doug Maguire too */
    if( plant_ptr->tht <= 0.5)                                       
    {
        //*height_growth = 0.20;
        plant_ptr->tht_growth = 0.20;
        return;
    }
    

    /* set the coefficients for the model  base b0-b5, veg v1 v2, tpa d1 d2, relht h0 and h1 */
    /****************************************************************************************/
    //if ( c_ptr->plant_type == CONIFER || c_ptr->plant_type == HARDWOOD ) 
    if ( coeff_ptr->type == CONIFER || coeff_ptr->type == HARDWOOD ) 
    {

        /* these are the cofficients that were approved for the first version of the model */
        //b1 = coeffs_ptr[0]; // = 11.75218;
        //b2 = coeffs_ptr[1]; // = -1.92283;
        //b3 = coeffs_ptr[2]; // = -0.0079;
        //b4 = coeffs_ptr[3]; // = -1.85283;
        //b5 = coeffs_ptr[4]; // = 1.019159;

        /* these are the cofficients that were approved for the first version of the model */
        //b1 = coeffs_ptr[0]; // = 13.2604;
        //b2 = coeffs_ptr[1]; // = -1.5678;
        //b3 = coeffs_ptr[2]; // = -0.0101;
        //b4 = coeffs_ptr[3]; // = -1.4203;
        //b5 = coeffs_ptr[4]; // = 0.1756;

        b1 = coeff_ptr->ht_growth[0]; // = 13.2604;
        b2 = coeff_ptr->ht_growth[1]; // = -1.5678;
        b3 = coeff_ptr->ht_growth[2]; // = -0.0101;
        b4 = coeff_ptr->ht_growth[3]; // = -1.4203;
        b5 = coeff_ptr->ht_growth[4]; // = 0.1756;


        /* compute psi first, which is used to compute the  */
        /* potential height growth                          */
        // STOCK_HT0 is the height of the plented stock?
		// x 	= Number of years since planting. Make X=30-stock age   */
        // since the stock age is 2.0?, make x=30 - stock age
        // stock_age = STOCK_AGE
        // x = 30.0 - STOCK_AGE, which in this case is assumed to be 2.0
        //get_psi( 28.0, FLEWELLING_ET_AL_2002_STOCK_HT0, flew_site, &psi );
        get_psi( 28.0, FLEWELLING_ET_AL_2002_STOCK_HT0, plot_ptr->site_30, &psi );


        /*Get growth effective age and compute 1-year h40 growth (in feet) this year*/
        get_delta_h_pot( psi, FLEWELLING_ET_AL_2002_STOCK_HT0, h40, &delta_h_pot );

        /* Convert to centimeters */
        delta_h_pot *= FT2CM;

        /*
  **************************************************************
  these are from the spreadsheet - CIPS PAI equations new.xls 
  **************************************************************

?H = 13.2604 + I*[ ?Hpot *(VegmodH)*(RHmodH)  ]                 
                
    I       =   0 if first growing season after planting; 1 otherwise
    ?Hpot       =   Potential top height growth from Flewelling et al. (2002) under weed-free conditions
    VegmodH     =   exp[–exp(-1.5678 – 0.0101*H)*Vegcov0.5]
    RHmodH      =   [ 1 – exp(-1.4203*Hrel0.1756) ]

    where delta H is height growth in cm            
    where delta Hpot is potential height growth in cm           
    where H is height in cm at beginning of growing season          
    where vegcov is summed weed cover in % i.e. 20% weed cover          
    where Hrel is relative height (ht of subject tree dividedby height of 40 largest per acre)          
    where I=0 if first growing season after planting; otherwise=1           

    dH          dHpot   ht  vegcov  Hrel            dH (in feet)
    62.93943054 115     50  20        1           2.064941947
    84.3091373  115    150  20        1           2.766047812
==> 99.05186605 115    400  20        1           3.249733138 <==
    41.5596345  115     50  80        1           1.363505069
    71.14222625 115    150  80        1           2.334062541
    97.6556062  115    400  80        1           3.203924088
*/

        /* compute the relative height modifier (eq. 3, page 32) */
        rel_height = (plant_ptr->tht * FT2CM ) / ( h40 * FT2CM );

        /* assign the temp variables for testing */
        //base_model      = 115.0; 
        //total_height    = 150.0 * CM2FT;
        //pct_veg_cov     = 80.0;
        //rel_height      = 1.0;
                
        /* compute the competeting vegetation modifier */
        //cv_mod_hg = exp( -exp( -1.5678 - 0.0101 * ( total_height * FT2CM ) ) * pow( pct_veg_cov, 0.5 ) );
        //cv_mod_hg = exp( -exp( b2 + b3 * ( total_height * FT2CM ) ) * pow( pct_veg_cov, 0.5 ) );
        cv_mod_hg = exp( -exp( b2 + b3 * ( plant_ptr->tht * FT2CM ) ) * pow( plot_ptr->shrub_pct_cover, 0.5 ) );

        /* compute the relative height modifier */
        //rh_mod_hg = 1.0 - exp( -1.4203 * pow( rel_height, 0.1756 ) );
        rh_mod_hg = 1.0 - exp( b4 * pow( rel_height, b5 ) );

        /* compute the final height growth estimate */
        //*height_growth = 13.2604 + after_first_season * ( delta_h_pot * cv_mod_hg * rh_mod_hg );
        //*height_growth = b1 + after_first_season * ( delta_h_pot * cv_mod_hg * rh_mod_hg );
        plant_ptr->tht_growth = b1 + after_first_season * ( delta_h_pot * cv_mod_hg * rh_mod_hg );
        
        /* and convert back into feet */
        //*height_growth *= CM2FT;
        plant_ptr->tht_growth *= CM2FT;
    }
    else if (coeff_ptr->type == SHRUB )   
    {
        //*height_growth = 0.0;
        plant_ptr->tht_growth = 0.0;
    }
    else /* its non stocked or forb */
    {
        //*height_growth = 0.0;
        plant_ptr->tht_growth = 0.0;
    }

    /* if predicted height is <0.5 then        */
//    if( total_height + *height_growth <= 0.5 )                       
//    {
//        *height_growth = (-1.0)*(total_height-0.51 ) ;
//    }

    if( plant_ptr->tht + plant_ptr->tht_growth <= 0.5 )                       
    {
        //*height_growth = (-1.0)*(total_height-0.51 ) ;
        plant_ptr->tht_growth = (-1.0f)*(plant_ptr->tht - 0.51f ) ;
    }



#ifdef _DEBUG
    /* assign the debugging/test output variables to the spares */
    /* I need a mapping between the spare and the debugging variable */
    /* lables for integer variables */
    plant_ptr->DBL_SPARE[0]         =   rel_height;
    plant_ptr->DBL_SPARE[1]         =   rh_mod_hg;
    plant_ptr->DBL_SPARE[2]         =   cv_mod_hg;
    plant_ptr->DBL_SPARE[3]         =   delta_h_pot;
    plant_ptr->DBL_SPARE[4]         =   psi;

    /* assign the integer spares and debugging labels */
    plant_ptr->INT_SPARE[0]         =   after_first_season;


#endif


    return;
}




/********************************************************************************/
/*                  cips_calc_dbh_growth                                        */
/********************************************************************************/
/*  Description :   cips_calc_dbh_growth                                        */   
/*  Author      :   Jeff D. Hamann & Douglas A. Maguire                         */
/*  Date        :   June 20, 2011                                               */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double         total_height     - total tree height                         */
/*  double         h40              - height of the 40 tallest stems            */
/*  double         current_dbh      - initial dbh (inches)                      */
/*  double         *pred_dbh_growth - predicted dbh growth (inches)             */
/*  vector<double> *coeffs_ptr      -   pointer to a vector of doubles that     */
/*                                      contain the coefficients for the        */
/*                                      functional species code                 */
/*  unsigned long  plant_type       - growth form of the plant record           */
/********************************************************************************/
/*  Formula : equation 4, page 32                                               */
/*  Source  : Center for Intensive Planted-Forest Silviculture,                 */
/*              2010 Annual Report                                              */
/*              http://www.fsl.orst.edu/cips/                                   */
/*  Coeffs  : DG                                                                */
/********************************************************************************/

static void cips_calc_dbh_growth(   
    unsigned long           *return_code,   
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  h40 )
{

    double  b1          = 0.0;
    double  b2          = 0.0;
    double  b3          = 0.0;
    double  b4          = 0.0;
    double  b5          = 0.0;
    double  b6          = 0.0;
    double  b7          = 0.0;
    double  b8          = 0.0;
    double  b9          = 0.0;
    double  rel_height  = 0.0; /* relative height = total_height / ht40*/
    double  rh_mod_dg   = 0.0;  /* equation 6, page 32 */
    double  cv_mod_dg   = 0.0;  /* equation 5, page 32 */
//    double  error_d     = 0.0;
//    double  site_model  = 0.0;
    double  base_model  = 0.0;

    *return_code=CONIFERS_SUCCESS;

    if( coeffs_ptr == NULL )
    {
        plant_ptr->dbh_growth    =   0.0f;
        *return_code        = INVALID_COEFF;
        return;
    }

    if( plant_ptr->tht <= 4.5f )
    {
        plant_ptr->dbh_growth    =   0.0f;
        return; 
    }

    b1  = coeffs_ptr->dbh_growth[0]; // b1	= 0.00000162;
    b2  = coeffs_ptr->dbh_growth[1]; // b2  = 0.097185;
    b3  = coeffs_ptr->dbh_growth[2]; // b3	= -0.00127;
    b4  = coeffs_ptr->dbh_growth[3]; // b4	= -3.43137;
    b5  = coeffs_ptr->dbh_growth[4]; // b5	= -0.00069;
    b6  = coeffs_ptr->dbh_growth[5]; // b6	= -4.84489;
    b7  = coeffs_ptr->dbh_growth[6]; // b7	= 0.003347;
    b8  = coeffs_ptr->dbh_growth[7]; // b8	= -4.52113;
    b9  = coeffs_ptr->dbh_growth[8]; // b9	= 4.644605;

    rel_height = (plant_ptr->tht * FT2CM ) / ( h40 * FT2CM );

    /* part 1 of the equation */    
    base_model = b1 * pow( plant_ptr->dbh * IN2MM, b2 ) * 
            exp( b3 * ( plant_ptr->dbh * IN2MM ) );

    /* veg cover modifier for the diameter growth model (veg_mod_dg) */
    /* check with mainwaring 3 2014*/
    cv_mod_dg       =   exp( -exp( b4 + b5 * ( plant_ptr->tht * FT2CM ) ) * 
        pow( plot_ptr->shrub_pct_cover, 0.5f ) );

    /* relative height growth modifier for the  */
    /* diameter growth model (rh_mod_dg)        */
    rh_mod_dg       =   exp( -exp( b6 + b7 * ( plant_ptr->tht * FT2CM ) ) * 
        pow( rel_height, b8 ) );

    /* then compute the predicted dbh growth for the plant  */
    /* updated equation provided by doug mainwaring         */
    /* August 21, 2011                                      */
    /* not right yet at 3.2014 addd basal area*/
    plant_ptr->dbh_growth =  base_model * cv_mod_dg  * rh_mod_dg  * 
        pow( ( plot_ptr->site_30 * FT2M ), b9 );  

    /* convert back into imperial units */
    plant_ptr->dbh_growth *= MM2IN;

    /* set the error flag to success, and return */
    *return_code            = CONIFERS_SUCCESS;

    if( plant_ptr->dbh_growth < 0.0f )
    {
        plant_ptr->dbh_growth    = 0.0f;
        *return_code        = CONIFERS_ERROR;
    }


#ifdef _DEBUG
/* assign the debugging/test output variables to the spares */

    plant_ptr->DBL_SPARE[5]         =   base_model;
    plant_ptr->DBL_SPARE[6]         =   cv_mod_dg;
    plant_ptr->DBL_SPARE[7]         =   rh_mod_dg;

#endif



}





/********************************************************************************/
/*                  calc_dbh_from_ht_and_veg_cov                                */
/********************************************************************************/
/*  Description :                                                               */
/* this function computes the dbh for the tree record, given the total height   */
/* and the proportion of the area covered by competing vegetation on the plot   */
/*  Author      :   Jeff D. Hamann and Douglas A. Maguire                       */
/*  Date        :   June 20, 2011                                               */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  pct_veg_cover               -   total height of the subject tree            */
/*  *pred_dbh                    -   predicted breast height diameter           */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : Not in Annual Report Document                                     */
/*  Source  : Center for Intensive Planted-Forest Silviculture,                 */
/*              2010 Annual Report                                              */
/*              http://www.fsl.orst.edu/cips/                                   */
/*  Coeffs  : DH                                                                */
/********************************************************************************/
//void cips_calc_dbh_from_ht_and_veg_cov( 
//    unsigned long *return_code,
//    double      total_height,
//    double      pct_veg_cover,
//    double      *pred_dbh,
//    double      *coeffs_ptr    )    /* dbh_ht_veg_cov */

static void cips_calc_dbh_from_ht_and_veg_cov(   
    unsigned long           *return_code,   
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr )


{

    double  b0  =   0.0;
    double  b1  =   0.0;
    double  b2  =   0.0;
    double  mse =   0.0;

    /* perform error check for the correct number of coeffs */
    if( coeffs_ptr == NULL )
    {
        *return_code = INVALID_COEFF;
        return;
    }

    if( plant_ptr->tht <= 4.5f )
    {
        *return_code = INVALID_INPUT_VAL;
        return;
    }

    /* dummy in the real coeffs for now */
    b0  =   coeffs_ptr->dbh_ht_veg_cov[0];      /* = -10.4003   */
    b1  =   coeffs_ptr->dbh_ht_veg_cov[1];      /* = 2.44412    */
    b2  =   coeffs_ptr->dbh_ht_veg_cov[2];      /* = 0.04612    */
    mse =   coeffs_ptr->dbh_ht_veg_cov[3];      /* = 0.061653   */

    /* this is a simple check to ensure the model can handle plots  */
    /* with no vegetation (debugging purposes)                      */
    if( plot_ptr->shrub_pct_cover <= 0.0f )
    {
        plot_ptr->shrub_pct_cover = 0.001f;
    }

    plant_ptr->dbh = exp( b0 + b1 * log( ( plant_ptr->tht * FT2CM ) ) + 
                            b2 * log( plot_ptr->shrub_pct_cover ) + 0.5 * mse );
    *return_code = CONIFERS_SUCCESS;

    plant_ptr->dbh *= MM2IN;

    if( plant_ptr->dbh < 0.0f )
    {
        plant_ptr->dbh  = 0.0f;
        *return_code    = CONIFERS_ERROR;
    }

}




/********************************************************************************/
/*                  calc_d6_from_ht_and_veg_cov                                 */
/********************************************************************************/
/*  Description :                                                               */
/* this function computes the basal diameter for the tree record, given the     */
/* total height and the proportion of the area covered by competing vegetation  */
/* on the plot                                                                  */
/*  Author      :   Jeff D. Hamann and Douglas A. Maguire                       */
/*  Date        :   June 20, 2011                                               */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  pct_veg_cover               -   total height of the subject tree            */
/*  *pred_dbh                    -   predicted breast height diameter           */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : Not in Annual Report Document                                     */
/*  Source  : Center for Intensive Planted-Forest Silviculture,                 */
/*              2010 Annual Report                                              */
/*              http://www.fsl.orst.edu/cips/                                   */
/*  Coeffs  : DH                                                                */
/********************************************************************************/
//void cips_calc_d6_from_ht_and_veg_cov( 
//    unsigned long           *return_code,
//    struct PLOT_RECORD      *plot_ptr,
//    struct PLANT_RECORD     *plant_ptr,
//    struct COEFFS_RECORD    *coeffs_ptr )
//{
//
//    double  part_a  = 0.0;
//    double  part_b  = 0.0;
//    double  X       = 0.0;

//    double  b1      = 0.0;
//    double  b2      = 0.0;
//    double  b3      = 0.0;
//    double  b4      = 0.0;
//    double  b5      = 0.0;
//    double  b6      = 0.0;
//    double  b7      = 0.0;
    

    /* perform error check for the correct number of coeffs */
//    if( coeffs_ptr == NULL )
//    {
//        *return_code = INVALID_COEFF;
//        return;
//    }

    //if( total_height * FT2CM <= 15.0 )
//    if( plant_ptr->tht * FT2CM <= 15.0 )
//    {
//        *return_code = INVALID_INPUT_VAL;
//        return;
//    }

//    b1  = coeffs_ptr->d6_ht_veg_cov[0];  /* = -3.4952 */
//    b2  = coeffs_ptr->d6_ht_veg_cov[1];  /* = 0.6412 */
//    b3  = coeffs_ptr->d6_ht_veg_cov[2];  /* = 0.2687 */
//    b4  = coeffs_ptr->d6_ht_veg_cov[3];  /* = 2.6949 */
//    b5  = coeffs_ptr->d6_ht_veg_cov[4];  /* = 0.2507 */
//    b6  = coeffs_ptr->d6_ht_veg_cov[5];  /* = -0.0370 */
//    b7  = coeffs_ptr->d6_ht_veg_cov[6];  /* = 0.970 */

    /* to prevent raising zero to a power */
    /* add a little tiny value to the pct_veg_cover if it's zero */
//    if( plot_ptr->shrub_pct_cover <= 0.000001f )
//    {
//        plot_ptr->shrub_pct_cover = 0.000001f;
//    }

    /* what the coeffs that need to go in here? */
//    X = 1.0 - pow( 15.0 / ( plant_ptr->tht * FT2CM ), 0.5 );
 
//    part_b = ( 15.0 / ( plant_ptr->tht * FT2CM ) );

    /* there's a bug in the equation provided by doug maguire */
//    part_a =    b1 * part_b + 
//                b2 * pow( part_b, 2.0f ) + 
//                b3 * log( part_b + 0.0001f) + 
//                b4 * pow( part_b, 0.5f );

    /* entered from contract */
//    plant_ptr->d6 = b5 * pow( plot_ptr->shrub_pct_cover, b6 ) * 
//        pow( ( plant_ptr->tht * FT2CM ), b7 ) * pow( X, part_a );
//    plant_ptr->d6 *= MM2IN;
//    *return_code = CONIFERS_SUCCESS;

//    if( plant_ptr->d6 < 0.0f )
//    {
//        plant_ptr->d6   = 0.0f;
//        *return_code    = CONIFERS_ERROR;
//    }


//}




/********************************************************************************/
/*                  calc_d12_from_ht_and_veg_cov                                */
/********************************************************************************/
/*  Description :                                                               */
/* this function computes the diameter at 30 cm for the tree record, given the  */
/* total height and the proportion of the area covered by competing vegetation  */
/* on the plot                                                                  */
/*  Author      :   Jeff D. Hamann and Douglas A. Maguire                       */
/*  Date        :   June 20, 2011                                               */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  pct_veg_cover               -   total height of the subject tree            */
/*  *pred_d12                   -   predicted d30 diameter                      */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : Not in Annual Report Document                                     */
/*  Source  : Center for Intensive Planted-Forest Silviculture,                 */
/*              2010 Annual Report                                              */
/*              http://www.fsl.orst.edu/cips/                                   */
/*  Coeffs  : DH                                                                */
/********************************************************************************/

//void cips_calc_d12_from_ht_and_veg_cov(
//    unsigned long           *return_code,
//    struct PLOT_RECORD      *plot_ptr,
//    struct PLANT_RECORD     *plant_ptr,
//    struct COEFFS_RECORD    *coeffs_ptr )
//{

//    double  part_a  = 0.0;
//    double  X       = 0.0;


//    double  b1      = 0.0;
//    double  b2      = 0.0;
//    double  b3      = 0.0;
//    double  b4      = 0.0;
//    double  b5      = 0.0;
//    double  b6      = 0.0;
//    double  b7      = 0.0;

    /* perform error check for the correct number of coeffs */
//    if( coeffs_ptr == NULL )
//    {
//        *return_code = INVALID_COEFF;
//        return;
//    }

    /* does this really apply? */
    /* question: what are the ranges on the inputs this function */
//    if( plant_ptr->tht * FT2CM <= 30.0f )
//    {
//        *return_code = INVALID_INPUT_VAL;
//        return;
//    }


//    b1  = coeffs_ptr->d12_ht_veg_cov[0];  /* = -3.4952 */
//    b2  = coeffs_ptr->d12_ht_veg_cov[1];  /* = 0.6412 */
//    b3  = coeffs_ptr->d12_ht_veg_cov[2];  /* = 0.2687 */
//    b4  = coeffs_ptr->d12_ht_veg_cov[3];  /* = 2.6949 */
//    b5  = coeffs_ptr->d12_ht_veg_cov[4];  /* = 0.2507 */
//    b6  = coeffs_ptr->d12_ht_veg_cov[5];  /* = -0.0370 */
//    b7  = coeffs_ptr->d12_ht_veg_cov[6];  /* = 0.970 */

    /* to prevent raising zero to a power */
    /* add a little tiny value to the pct_veg_cover if it's zero */
//    if( plot_ptr->shrub_pct_cover <= 0.000001f )
//    {
//        plot_ptr->shrub_pct_cover = 0.000001f;
//    }

    /* what the coeffs that need to go in here? */
//    X = 1.0f - pow( 30.0f / ( plant_ptr->tht * FT2CM ), 0.5f );
 
    /* there's a bug in the equation provided by doug maguire */
//    part_a = b1 * ( 30.0f / ( plant_ptr->tht * FT2CM ) ) + 
//                b2 * pow( 30.0f / ( plant_ptr->tht * FT2CM ), 2.0f ) + 
//                b3 * log( ( 30.0f / ( plant_ptr->tht * FT2CM ) ) + 0.0001f ) + 
//                b4 * pow( 30.0f / ( plant_ptr->tht * FT2CM ), 0.5f );

    /* entered from contract */
    //*pred_d12 = 0.2507 * pow( pct_veg_cover, -0.0370 ) * pow( ( total_height * FT2CM ), 0.970 ) * pow( X, part_a );
    //*pred_d12 = b5 * pow( pct_veg_cover, b6 ) * pow( ( total_height * FT2CM ), b7 ) * pow( X, part_a );
    //*return_code = CONIFERS_SUCCESS;

//    plant_ptr->d12 = b5 * pow( plot_ptr->shrub_pct_cover, b6 ) * 
//        pow( ( plant_ptr->tht * FT2CM ), b7 ) * pow( X, part_a );

//    plant_ptr->d12 *= MM2IN;
    
//    *return_code = CONIFERS_SUCCESS;

//    if( plant_ptr->d12 < 0.0f )
//    {
//        plant_ptr->d12       = 0.0f;
//        *return_code    = CONIFERS_ERROR;
//    }
//
//
//}



/********************************************************************************/
/*                  cips_calc_crown_ratio                                       */
/********************************************************************************/
/*  Description :   calc_crown_ratio                                            */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   May 9, 2008                                                 */
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
/*  Formula : CL =b0*H**b1 *exp( b3 d6/height)                                  */ 
/*  Source  : Ritchie May 2008                                                  */
/*  Coeffs  : CR                                                                */
/********************************************************************************/
static void cips_calc_crown_ratio(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr )
{

    double  b0;
    double  b1;
    double  newCR;           

    *return_code = CONIFERS_SUCCESS;
    if( coeffs_ptr == NULL )
    {
        *return_code = INVALID_COEFF;
        plant_ptr->cr = 0.0;
        return;
    }

    if( plant_ptr->tht <= 0.0f )
    {
        *return_code = INVALID_INPUT_VAL;
        plant_ptr->cr = 0.0f;
        return;
    }

    b0 = coeffs_ptr->crown_ratio[0]; /* -0.00079 */
    b1 = coeffs_ptr->crown_ratio[1]; /* 0.8935   */
    
    newCR = 1.0 - ( 1.0 - exp( b0 * pow( plant_ptr->tht * FT2CM, b1 ) ) );

    if( newCR < 0.0 )
    {
        newCR = 0.0;
    }
    if( newCR > 1.0 )
    {
        newCR = 1.0;
    }

    plant_ptr->cr = newCR;

    /* no conversion needed */

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
static void cips_calc_max_crown_width(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr )
{
    double  b0;
    double  b1;
    double  b2;

    if( coeffs_ptr == NULL )
    {
        plant_ptr->max_crown_width  = 0.0f;
        *return_code                = INVALID_COEFF;
        return;
    }
    
    if( plant_ptr->dbh < 0.0f )
    {
        *return_code                = INVALID_INPUT_VAL;
        plant_ptr->max_crown_width  = 0.0;
        return;
    }

    b0 = coeffs_ptr->max_crown_width[0];
    b1 = coeffs_ptr->max_crown_width[1];
    b2 = coeffs_ptr->max_crown_width[2];
    
    if( plant_ptr->tht <= 4.5f )
    {
        plant_ptr->max_crown_width = b0 * (plant_ptr->tht / 4.5f);
    }
    else
    {
        plant_ptr->max_crown_width = b0 + b1 * plant_ptr->dbh + b2 * 
            plant_ptr->dbh * plant_ptr->dbh;
    }

    *return_code = CONIFERS_SUCCESS;

    if( plant_ptr->max_crown_width < 0.0f )
    {
        plant_ptr->max_crown_width   = 0.0f;
        *return_code            = CONIFERS_ERROR;
    }


}



/********************************************************************************/
/*                  smc_calc_crown_width     S1                                 */
/********************************************************************************/
/*  Description :   smc_calc_crown_width                                        */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   November 30, 2007                                           */
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
static void cips_calc_crown_width(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr )
{

    double  b0;
    double  b1;
    double  b2;

        /* check for valid height */
        if(plant_ptr->tht <= 0.0f)
        {
            *return_code            = INVALID_INPUT_VAL;
            plant_ptr->crown_width = 0.0f;
            plant_ptr->crown_area  = 0.0f;
            return;
        }

        /* check for valid coefficients */
        if( coeffs_ptr == NULL )
        {
            *return_code        = INVALID_COEFF;
            plant_ptr->crown_width = 0.0f;
            plant_ptr->crown_area  = 0.0f;
            return;
        }

        b0 = coeffs_ptr->crown_width[0];
        b1 = coeffs_ptr->crown_width[1];
        b2 = coeffs_ptr->crown_width[2];

        if( plant_ptr->tht < 0.51f )
        {
            plant_ptr->crown_width = 0.25f;
            plant_ptr->crown_area  = 0.04908739f;
        }
        else
        {
            plant_ptr->crown_area = exp( b0 + 
                                    b1 * log( plant_ptr->d6_area * 144.0f ) + 
                                    b2 * log( plant_ptr->tht ) );

            if(coeffs_ptr->type == SHRUB)
			{
				if( plant_ptr->crown_area > 50.1f )
				{
					plant_ptr->crown_area = 50.1f ; /* limit crown width to 8 feet for shrubs */
				}
			}
			else
			{
				if( plant_ptr->crown_area > 2827.0f )
				{
			        plant_ptr->crown_area = 2827.0f ; /* limit crown width to 60 feet */
				}
			}

			plant_ptr->crown_width = sqrt( plant_ptr->crown_area * ONE_OVER_PI * 4.0f); 
        }

        *return_code = CONIFERS_SUCCESS;

        if( plant_ptr->crown_width < 0.0f )
        {
            plant_ptr->crown_width   = 0.0f;
            plant_ptr->crown_area    = 0.0f;
            *return_code        = CONIFERS_ERROR;
			return;
        }
}


/********************************************************************************/
/*                  cips_calc_endemic_mortality                                 */
/********************************************************************************/
/*  Description :   cips_calc_endemic_mortality                                 */   
/*  Author      :   Jeff D. Hamann & Douglas A. Maguire                         */
/*  Date        :   June 20, 2011                                               */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double         expansion_factor - current expansion factor                  */
/*  double         total_height     - total tree height                         */
/*  double         h40              - height of the 40 tallest stems            */
/*  double         current_dbh      - initial dbh (inches)                      */
/*  double         *pred_mortality  - predicted mortality                       */
/*  vector<double> *coeffs_ptr      -   pointer to a vector of doubles that     */
/*                                      contain the coefficients for the        */
/*                                      functional species code                 */
/*  unsigned long  plant_type       - growth form of the plant record           */
/********************************************************************************/
/*  Formula : equation 2, page 38                                               */
/*  Source  : Center for Intensive Planted-Forest Silviculture,                 */
/*              2010 Annual Report                                              */
/*              http://www.fsl.orst.edu/cips/                                   */
/*  Coeffs  : Mortality                                                         */
/********************************************************************************/
static void cips_calc_endemic_mortality(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  h40 )
{

    double  b0;
    double  b1;
    double  b2;
    double  b3;
    double  prob_of_mortality;
    double  Y;
    double  rel_height;

    /*  initialize parameters    */
    plant_ptr->expf_change  = 0.0f;

    if( coeffs_ptr == NULL )
    {
        plant_ptr->expf_change  = 0.0f;
        *return_code            = INVALID_COEFF;
        return;
    }
    
    b0  = coeffs_ptr->mortality[0];  // = -1.3409
    b1  = coeffs_ptr->mortality[1];  // = -0.00328
    b2  = coeffs_ptr->mortality[2];  // = 5.0445
    b3  = coeffs_ptr->mortality[3];  // = 0.1634

    if( coeffs_ptr->type == CONIFER )
    {
        rel_height = ( plant_ptr->tht * FT2CM ) / ( h40 * FT2CM );

        /*   this will calculate the number of trees to die using the endemic mortality and current mort if any */
        /* todo: add switch for when pct_veg_cover == 0.0, log( 0.0 ) = -inf */
        if( plot_ptr->shrub_pct_cover <= 0.0f )
        {
            plot_ptr->shrub_pct_cover = 0.001f;
        }

        /* used for testing */
        //total_height    = 175.0 * CM2FT;
        //rel_height      = 1.0;
        //pct_veg_cover   = 50.0;
        Y = b0 + b1 * ( plant_ptr->tht * FT2CM ) + b2 * rel_height + b3 * log( plot_ptr->shrub_pct_cover );

        prob_of_mortality = exp( Y ) / ( 1.0 + exp( Y ) );    

        plant_ptr->expf_change = plant_ptr->expf * prob_of_mortality;
    }
    else
    {
        plant_ptr->expf_change     = 0.0f;
    }
    
    *return_code        = CONIFERS_SUCCESS;

    if( plant_ptr->expf_change < 0.0f )
    {
        plant_ptr->expf_change     = 0.0f;
        *return_code        = CONIFERS_ERROR;
    }


#ifdef _DEBUG
    /* assign the debugging/test output variables to the spares */
    plant_ptr->DBL_SPARE[8]         =   prob_of_mortality;

#endif



}



/********************************************************************************/
/*                  cips_calc_cw_growth                                         */
/********************************************************************************/
/* these values are plant level variables                                       */
static void cips_calc_cw_growth(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr,
    unsigned long           plantation_age )
{

    double  initial_pct_cov;
    double  temp_cover_growth;

	double  temp_cover;
	double  temp_cwg;
	double  c0;
	double  c1;
	double  c2;
	double  c3;
	double  c4;
	double  c5;
	double  c6;
	//double  c7;  // unused jan 2014 mwr;
	//double  c8;  // unused jan 2014 mwr;
	//double  c9;  // unused jan 2014 mwr;
	//double  c10; // unused jan 2014 mwr;

    double  ca_conifers;
	double  ca_hardwoods;
	double  ca_shrubs;

	*return_code	= CONIFERS_SUCCESS;

    if( coeffs_ptr == NULL )
    {
        plant_ptr->cw_growth    = 0.0;
#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
        *return_code        = INVALID_COEFF;
        return;
    }

    if( plot_ptr->ca_c < 0.0f || 
        plot_ptr->ca_h < 0.0f || 
        plot_ptr->ca_s < 0.0f )
    {
 	    plant_ptr->cw_growth	=   0.0f;
#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
	    *return_code	        =   INVALID_INPUT_VAL;
	    return;
    }

    if( plant_ptr->crown_width <= 0.0f )
    {
 	    plant_ptr->cw_growth	=   0.0f;
#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
	    *return_code	        =   INVALID_INPUT_VAL;
	    return;
    }

    if( plant_ptr->tht_growth < 0.0f ) 
    {
        plant_ptr->cw_growth    = 0.0f;
#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
	    *return_code	        = CONIFERS_SUCCESS;
	    return;
    }

    temp_cwg    = 0.0f;

    /* then use the smc version for trees       */
    /* this is originally from the SMC model    */
    if( coeffs_ptr->type == CONIFER || coeffs_ptr->type == HARDWOOD )
    {

        ca_conifers     = ( plot_ptr->ca_c / SQ_FT_PER_ACRE );
        ca_hardwoods    = ( plot_ptr->ca_h / SQ_FT_PER_ACRE );
        ca_shrubs       = ( plot_ptr->ca_s / SQ_FT_PER_ACRE );

        c0  = coeffs_ptr->cw_growth[0];
        c1  = coeffs_ptr->cw_growth[1];
        c2  = coeffs_ptr->cw_growth[2];
        c3	= coeffs_ptr->cw_growth[3];
        c4	= coeffs_ptr->cw_growth[4];
        c5	= coeffs_ptr->cw_growth[5];
        c6	= coeffs_ptr->cw_growth[6];
        //c7	= coeffs_ptr->cw_growth[7];   // unused jan 2014 mwr;
        //c8	= coeffs_ptr->cw_growth[8];   // unused jan 2014 mwr;
        //c9  = coeffs_ptr->cw_growth[9];   // unused jan 2014 mwr;
        //c10 = coeffs_ptr->cw_growth[10];  // unused jan 2014 mwr;

        /* this is the old cw growth model from the smc sim */
///  This needs to be checked mwr Jan 2014 pending check;
        temp_cwg = pow( plant_ptr->tht_growth, c1 )
                        * ( c0 
                        + c2 * sqrt( plant_ptr->crown_width ) 
                        + c3 * ca_conifers * ca_conifers 
                        + c4 * ca_hardwoods * ca_hardwoods 
                        + c5 * ca_shrubs * ca_shrubs 
                        + c6 * log( plant_ptr->crown_width ) );

      if( temp_cwg < 0.0f ) /* if pred growth is negative */
	  {
		 //*pred_cw_growth    = 0.10;  
          plant_ptr->cw_growth  = 0.10f;    /* Growth will be set to 0.10 */
#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
		 *return_code           = CONIFERS_ERROR;
		 return;
	  }

        plant_ptr->cw_growth    = temp_cwg;
#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
	    *return_code              = CONIFERS_SUCCESS;
	    return;
    }
    else if ( coeffs_ptr->type == SHRUB ) 
    {

        /* compute the converted vegetation cover (cv) where the */
        /* vc should be equal to ca_shrubs if one plant per plot */
        //initial_pct_cov = 100.0 * (MY_PI/(4.0*SQ_FT_PER_ACRE))*(crown_width*crown_width)*expf;
        initial_pct_cov = 100.0f * (MY_PI/( 4.0f * SQ_FT_PER_ACRE ) ) * 
                                    (plant_ptr->crown_width * plant_ptr->crown_width) * 
                                    plant_ptr->expf;
      
        /* added from doug maiwaring's fits */
        /* these are the second half of the coeffs */
        c0 = coeffs_ptr->cw_growth[11]; // =   125.8;                      
        c1 = coeffs_ptr->cw_growth[12]; // = -1.1124;            
        c2 = coeffs_ptr->cw_growth[13]; // = 0.07498;                        
        c3 = coeffs_ptr->cw_growth[14]; // = -8.3251;            
        c4 = coeffs_ptr->cw_growth[15]; // = -0.06569;                           

        /* you need to conver the change in total cover to the change in crown width    */
        /* for the plant                                                                */
        /* in this equation conifer means douglas-fir                                   */
        /* these values are normally negative                                           */
        //temp_cover_growth = c0 + 
        //                    ( c1 + c2 * plantation_age ) * plot_ptr->shrub_pct_cover + 
        //                    ( c3 + c4 * ( conifer_d12_basal_area * FTAC2M2HA ) ) * plantation_age;
        

        /* THIS IS WHERE THE PROBLEM IS */
        /* use the plot level variable for the d12 basal area in conifers */
        temp_cover_growth = c0 + 
                            ( c1 + c2 * plantation_age ) * plot_ptr->shrub_pct_cover + 
                            ( c3 + c4 * ( plot_ptr->d12ba_c * FTAC2M2HA ) ) * plantation_age;

        /* meaning that the results of the equation predicted a 3.3459999999999894 percent */
        /* increase in the percent cover. This value needs to be translated into the change in crown width */
        /* which is the variable that's tracked through the simulator */

        /* crown width growth */
        //temp_cover = initial_pct_cov + temp_cover_growth;
        temp_cover = initial_pct_cov + temp_cover_growth;

        /* if there are plants, then compute the temporary crown width growth */
        if( plant_ptr->expf > 0.0f )
        {
            temp_cwg = sqrt((temp_cover/100.0f)*(SQ_FT_PER_ACRE*4.0f/(MY_PI*plant_ptr->expf))) - plant_ptr->crown_width;
        }

        /* if the crown width for the plant will grow to be */
        /* over 8 feet in diameter, then keep shrubs from   */
        /* getting over 8 feet in diameter                  */
        if( temp_cwg + plant_ptr->crown_width > 8.0f ) 
        {
            //plant_ptr->crown_width  = 0.0f;
            plant_ptr->cw_growth  = 0.0f;
#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
            *return_code            = CONIFERS_SUCCESS;
            return;
        }

        
        if( temp_cwg < 0.0 ) /* if cover is decreasing make a change in mortality */
        {
            //*pred_mortality = plant_ptr->expf - (temp_cover/100.0)*
            //                    (SQ_FT_PER_ACRE*4.0 / (MY_PI*plant_ptr->crown_width*plant_ptr->crown_width));

            //plant_ptr->cw_growth = temp_cwg;
            
            plant_ptr->expf_change = plant_ptr->expf - (temp_cover/100.0f)*
                                (SQ_FT_PER_ACRE*4.0f / (MY_PI*plant_ptr->crown_width*plant_ptr->crown_width));

            plant_ptr->cw_growth = temp_cwg;

#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
            *return_code    = CONIFERS_SUCCESS;
            return;
        }


	  plant_ptr->cw_growth  = temp_cwg;


#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif

        *return_code     = CONIFERS_SUCCESS;
	  return;

    }
	else  /* then you don't have a valid plant type for this variant */
	{
      plant_ptr->cw_growth   =   0.0f;
#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif
	  *return_code      =   INVALID_PLANT_TYPE;
	  return;
	}


#ifdef _DEBUG
    plant_ptr->DBL_SPARE[9]         =   plant_ptr->cw_growth;
#endif

}





/* this function computes the potential height growth, given the    */
/* average height of the planted seedlings, the current top height  */
/* and the psi value, which is computed elsewhere                   */
void get_delta_h_pot (
    double    psi,          /* PSI  = PSI computed elsewhere                       */
    double    ht0,          /* HT0  = average height of the planted stock (ft)     */
    double    htx,          /* HTx = current top height                            */
    double   *delta_output  /* output estimate of potential 1-year growth   */
)
{
    double gea      =   0.0f; /* temp variable for growth effective age */
    double htop1    =   0.0f; /* temp variable for top height           */ 

    /* compute the growth effective age */
    /* equation x, page p?.             */
    get_gea(psi, ht0, htx, &gea);
    
    /* Compute one year growth with current psi and growth effective age */
    /* Put estimate of next year's htop in last_htop */
    flewelling_site_index( psi, ht0, gea + 1.0f, &htop1 );

    /* Store the expected one-year growth in delta_output */
    *delta_output = htop1 - htx;
    return;
}


/* this function computes the growth effective age */
void get_gea (
    double    psi, /* PSI  = PSI computed elsewhere                       */
    double    ht0, /* HT0  = average height of the planted stock (ft)     */
    double    htx, /* HTx = current top height                            */
    double   *gea_output /* output growth effective age                 */
)
{
    double x = 50.0f, last_age = 60.0f; /*Will store age to compute growth effective age */
    double tmp_htop = 0.0f, last_htop = 0.0f; /*Will store height output from _flewsi*/ 
    int place = -1; /*Store current decimal place*/

    /* If supplied height is less than stock height, give stock age as gea*/
    if ( htx < ht0 ) 
    {
        *gea_output = 2.0f;
        return;
    }

    /*Compute growth effective age (from findgea SAS macro)*/
    flewelling_site_index( psi, ht0, x, &tmp_htop );

    /*Decrement by 10, then 1, then 0.1, then 0.01, etc*/
    while ( place <= 3 ) 
    {
        /*Drop psi until htop drops below our observed height*/
        while ( tmp_htop > htx ) 
        {
            last_age = x;
            last_htop = tmp_htop;
            INCDOWN(x, place);
            flewelling_site_index( psi, ht0, x, &tmp_htop );
        }

        x = last_age;
        tmp_htop = last_htop;
        place += 1;

    }

    *gea_output = x - 0.0005f;
    return;
}


/* x    = Number of years since planting. Make X=30-stock age   */
/* HT0  = average height of the planted stock (ft)      */
/* HTX  = top height at age x                           */
/* output estimate of psi                               */

void get_psi (
   double  x,           /* x    = Number of years since planting. Make X=30-stock age   */
   double ht0,          /* HT0  = average height of the planted stock (ft)      */
   double htx,          /* HTX  = top height at age x                           */
   double *psi_output   /* output estimate of psi                               */
)
{
    /*Start with a real high value of psi (5)*/
    double  cur_psi = 5.0f;
    double  last_psi = 0.0f;
    double  htop = 0.0;
    double  last_htop = 0.0;
    int     place = 0;

    /*If supplied height is less than stock height, output really low psi
     *which should drop growth to near 0, indicating a problem */
    if ( htx < ht0 ) 
    {
        *psi_output = 0.01f;
        return;
    }

    /* compute the initial value for htop */ 
    flewelling_site_index( cur_psi, ht0, x, &htop );
    
    /* solve htop to six decimal places             */
    /* Decrement by 1, then 0.1, then 0.01, etc     */
    while ( place <= 6 ) 
    {
        /*Drop psi until htop drops below our observed height*/
        while ( htop > htx ) 
        {
            last_psi = cur_psi;
            last_htop = htop;
            INCDOWN(cur_psi, place);
            flewelling_site_index(cur_psi, ht0, x, &htop);
        }
        
        cur_psi = last_psi;
        htop = last_htop;
        place += 1;
    }
    *psi_output = cur_psi - 0.0000005f;
    return;
}



/*********************************************************************************/
/*
this is the original source code from CIPS
* Compute site index (age in ft age total age=30 yrs) for a given PSI;
* The SI is actually defined by PSI, so use this to compute the normal SI value
	of height at 30 yrs total age;
* No local input variables for the macro, so embed this code into a sas dataset;
* All variables are years of total age, or ft;
* You need to supply
* 	PSI	= PSI computed elsewhere, 
	HT0 = average height of the planted stock (ft)
	x 	= Number of years since planting. Make X=30-stock age;

*  Compute site index (age in ft age total age=30 yrs) for a given PSI;
*   The SI is actually defined by PSI, so use this to compute the normal SI value*   of height at 30 yrs total age;
*   No local input variables for the macro, so embed this code into a sas dataset;
*   All variables are years of total age, or ft; 
*/
/*********************************************************************************/
void flewelling_site_index(
		double  psi,    /* PSI	= PSI computed elsewhere                                */
		double  ht0,    /* HT0  = average height of the planted stock (ft)              */
		double  x,      /* x    = Number of years since planting. Make X=30-stock age   */
        double  *htop_output
      )
{

    /*Variable names taken from original Fortran code in DF_PSET.FOR */
    double fp[14] = { 0, 6464.0f, -1691.0f, -29.23f, 7.510f, 0.9075f,
                     0.1788f, -102.4f, 31.93f, 39.67f, 66.58f, -42.79f,
                     16.57f, 0.3505f };
    double term, temp, z;
    double b1, c, xk1, xk2, yk1, yk2, xstraight;
    double top, shape;
    double lambda1, ll, alpha, beta, lambda2;
    double htop  = 0.0;

    /*Check for early age, return simply stock height at planting*/
    if ( x < 0.0 ) {
        *htop_output = FLEWELLING_ET_AL_2002_STOCK_HT0;
        return;
    }

    *htop_output = 0.0;

    /*Growth rate at origin (x=0, age=2) is b1 * PSI*/
    term = BOUND(fp[3] + fp[4] * psi, -8.0f, 8.0f);

    b1  = exp( term ) / ( 1 + exp( term ) );
    
    /*Inflection is at x= xk*/
    temp = fp[1] + fp[2] * psi;
    if ( temp > 0.0f )
    {
        temp = pow( temp, fp[13] );   
    }

    xk1 = MAX( 1.0f, temp ); 

    /*Shape of the below-inflection curve*/
    c  = fp[5] + fp[6] * psi;

    /*YK is value of y at inflection*/
    yk1 = ht0 + xk1*psi * (1.0f - (1.0f - b1) / (c + 1.0f) );

    /*insert the straight segment here*/
    xstraight = MAX( 0.0f , fp[7] + fp[8] * psi );

    xk2 = xk1 + xstraight;
    yk2 = yk1 + xstraight * psi;

    /*SLOPE  dY/dZ at Inflection  YPK
      if dZ/dX =1  then YPK = PSI
      else PSI = dY/dZ * dZ/dX
      (or YPK = PSI / (dZ/dX) ) = f11 + f12 * psi;*/

    /*TOP is the asymptote (ht at an infinite age).*/
    top = MAX( fp[9] * psi + fp[10], yk2 + 15.0f);

    shape =  BOUND(fp[11] + fp[12] * psi, -8.0f, 8.0f); 
    lambda1 = 0.05f + ( exp( 2.0f * psi / ( yk2 - top ) ) - 0.05f ) * 
              exp( shape ) / (1.0f + exp( shape ) );
    ll = log(lambda1);
    alpha = -psi * psi / 
            ( ll * ll * ( yk2 - top ) - 2.0f * ll * psi);
    beta = yk2 - top - alpha;
    lambda2 = exp( ( psi - alpha * ll ) / beta );

    /*Coefficients are computed as in DF_PSET.FOR, now use these to get htop*/
    if( x < xk1 )
    {
        htop = ht0 + psi * 
            ( x + ( 1 - b1 ) * xk1 / ( c + 1 ) * 
            ( pow( 1 - x / xk1 , ( c + 1 ) ) - 1 ) );
    }

    if( xk1 <= x && x <= xk2 )
    {
        htop = yk1 + ( x - xk1 ) * psi;
    }

    if( x > xk2 ) 
    {
    	z = x - xk2;
    	htop = yk2 + 
               alpha * ( ( pow( lambda1, z ) ) -1 ) + 
               beta *  ( ( pow( lambda2, z ) ) -1 );
    }


    /* assign the output variable */
    *htop_output = htop;
    return;
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
static void cips_calc_exp_from_cover_and_ca(
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr )
{
    double  cover_total=0.0;
    *return_code        = CONIFERS_ERROR;
    
    /* Cover must be between 0 and 100 per record   */
    if( plot_ptr->shrub_pct_cover < 0.0f || plot_ptr->shrub_pct_cover > 100.0f )
    {
        *return_code    = CONIFERS_ERROR;   
        plant_ptr->expf      = 0.0f;
        return;
    }

    /* crown area must be > 0 to avoid bad things   */
    if( plant_ptr->crown_area <= 0.0f )
    {
        *return_code    = CONIFERS_ERROR;   
        plant_ptr->expf      = 0.0f;
        return;
    }

    /* This is the total cover in sq ft from %cvr */
    cover_total     =   SQ_FT_PER_ACRE * plant_ptr->pct_cover * 0.01;
    plant_ptr->expf =   cover_total / plant_ptr->crown_area;           
    *return_code    =   CONIFERS_SUCCESS;
}



/********************************************************************************/
/* cips_impute                                                                   */
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
void cips_impute( 
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

    //double        cait_c;  // unused removed mwr jan 2014;
    //double        cait_h;  // unused removed mwr jan 2014;
    //double        cait_s;  // unused removed mwr jan 2014;

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
    /* first, impute the crown_area, needed for the simulator to    */
    /* compute the shrub_pct_cover (needed to compute d6, d12, and  */
    /* dbh?                                                         */
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        /* update the extra data for the tree                     */
        /* TODO: we should probably put this into another function */
        //plant_ptr->basal_area    = plant_ptr->dbh * plant_ptr->dbh * FC_I;
        //plant_ptr->d6_area       = plant_ptr->d6 * plant_ptr->d6 * FC_I;
        //plant_ptr->d12_area       = plant_ptr->d12 * plant_ptr->d12 * FC_I;
        
        /* you should only be tallying the shrubs */
        plant_ptr->crown_area    = plant_ptr->crown_width * 
                                    plant_ptr->crown_width * MY_PI / 4.0;
        plant_ptr->pct_cover     = 100.0 * plant_ptr->expf * 
                                    plant_ptr->crown_area/SQ_FT_PER_ACRE;
    }
    
    /* for this pass, we will only have the pct_shrub_cover */
    calc_plot_stats_2(  return_code, 
                        n_species,
                        species_ptr,
                        n_coeffs,
                        coeffs_ptr,
                        n_plants,
                        plants_ptr,
                        n_points,
                        plots_ptr );
    

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
	if( plant_ptr->tht >= 0.50f )
	{
	  if( is_tree( c_ptr ) )
	  {
		if( plant_ptr->d6 == 0.0f )
		{
		    if( plant_ptr->dbh > 0.0 && plant_ptr->tht >4.5f )
		    {
			    *return_code = CONIFERS_SUCCESS;
			
			  /* fill in the D6 */
			  /* set the pointer for the current plot that the plant is on */
			  plot_ptr = get_plot( plant_ptr->plot, n_points, plots_ptr );                  


                /* this function needs the veg cover to be correct      */
                /* so you better have figured out what it's supposed    */
                /* be BEFORE you get here.                              */
                //cips_calc_d6_from_ht_and_veg_cov( return_code,
                //                                    plot_ptr,
                //                                    plant_ptr,
                //                                    c_ptr );

                cips_calc_dob_hi( return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    15.0 * CM2FT,
                                    &plant_ptr->d6 );

			  /* if the total height of the plant > 30 cm         */
			  /* then fill in the missing d12 observation too     */
			  if( plant_ptr->tht > 30.0f * CM2FT )
				{


//                cips_calc_d12_from_ht_and_veg_cov( return_code,
//                                                    plot_ptr,
//                                                    plant_ptr,
//                                                    c_ptr );

                cips_calc_dob_hi( return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    30.0 * CM2FT,
                                    &plant_ptr->d12 );

				}
				
			    
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

			  
			      /* function not working */
			      /* set the pointer for the current plot that the plant is on */
			      plot_ptr = get_plot( plant_ptr->plot, n_points, plots_ptr );

//                  cips_calc_d6_from_ht_and_veg_cov( return_code,
//                                    plot_ptr,
//                                    plant_ptr,
//                                    c_ptr );

                cips_calc_dob_hi( return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    15.0 * CM2FT,
                                    &plant_ptr->d6 );
                                  
			    if( *return_code != CONIFERS_SUCCESS)
			      {
				plant_ptr->errors |= E_INVALID_D6;
				error_count += 1;
			      }
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
            
            plot_ptr = get_plot( plant_ptr->plot, n_points, plots_ptr );

//            cips_calc_d12_from_ht_and_veg_cov(  return_code,
//                                                plot_ptr,
//                                                plant_ptr,
//                                                c_ptr );

                cips_calc_dob_hi( return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    30.0 * CM2FT,
                                    &plant_ptr->d12 );
        
        
        }

          /* this is the code for computing the missing dbh when the tree */
          /* is taller than 4.5 feet                                        */
		if( plant_ptr->d6 > 0.0 && plant_ptr->dbh == 0.0 && plant_ptr->tht > 4.5)
		  {
		    *return_code = CONIFERS_SUCCESS;
		    
	      switch (variant)
			{
			  /* function not checked */
			  /* set the pointer for the current plot that the plant is on */
			  plot_ptr = get_plot( plant_ptr->plot, n_points, plots_ptr );
			    
			    cips_calc_dbh_from_ht_and_veg_cov(  return_code,
                                                    plot_ptr,
                                                    plant_ptr,
                                                    c_ptr );

//                cips_calc_dob_hi( return_code,
//                                    plot_ptr,
//                                    plant_ptr,
//                                    c_ptr,
//                                    137.0 * CM2FT,
//                                    &plant_ptr->dbh );
                
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
		      
			    /* set the pointer for the current plot that the plant is on */
			    plot_ptr = get_plot( plant_ptr->plot, n_points, plots_ptr );
			      
                /* this only applies to shrubs (is_shrub(c_ptr)==TRUE) */
//                cips_calc_d6_from_ht_and_veg_cov( return_code,
//                                                    plot_ptr,
//                                                    plant_ptr,
//                                                    c_ptr );

                cips_calc_dob_hi( return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr,
                                    15.0 * CM2FT,
                                    &plant_ptr->d6 );

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
		
		      cips_calc_crown_width(   return_code,
                                        plot_ptr,
                                        plant_ptr,
                                        c_ptr );

		    if( *return_code != CONIFERS_SUCCESS)
		      {
			plant_ptr->errors |= E_INVALID_CW;
			error_count += 1;
		      }
	      }

	    /* this should happen no matter what... */
	    plant_ptr->crown_area = plant_ptr->crown_width * plant_ptr->crown_width * MY_PI / 4.0;

	    /* fill in the expansion factors for those plant records    */
	    /* that don't have one filled in by using the sample        */
	    /* design records. Mostly, shrubs will have the expf filled */
	    /* by calling this function, but it should work for trees   */
	    if( plant_ptr->expf <= 0.0 )
	      {
              cips_calc_exp_from_cover_and_ca( return_code,
                                                plot_ptr,
                                                plant_ptr,
                                                c_ptr );

            if( *return_code != CONIFERS_SUCCESS)
		    {
		      plant_ptr->errors |= E_INVALID_EXPF;
		      error_count += 1;
		    }

	      }
    
    } /* END OF PASS ONE LOOP. */

    /* you need to recompute the plot stats */
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




    /* todo: this is where you need to put the udpated imputation functions */
    /* for d6 and d12 since they require plot_ptr->shrub_pct_cover.         */



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
            

		    //cait_c       =   cait[CONIFER];     // unused removed jan 2014 mwr;
		    //cait_h       =   cait[HARDWOOD];    // unused removed jan 2014 mwr;
		    //cait_s       =   cait[SHRUB];       // unused removed jan 2014 mwr;
		    *return_code = CONIFERS_SUCCESS;

			/* todo: make sure places where variants show up, that the functions are updated */
			/* not checked */
            cips_calc_crown_ratio(return_code,
                                    plot_ptr,
                                    plant_ptr,
                                    c_ptr );

        }                              
	  }

	if( *return_code != CONIFERS_SUCCESS)
	  {
	    plant_ptr->errors |= E_INVALID_CR;
	    error_count += 1;
	  }

	plant_ptr->max_crown_width = 0.0;

	  /* if it's a tree, compute the crown width */
	  if( is_tree( c_ptr ) )
	    {
	      *return_code = CONIFERS_SUCCESS;
     
		    cips_calc_max_crown_width(  return_code,
                                        plot_ptr,
                                        plant_ptr,
                                        c_ptr );

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



/* this is a diagnostic output file that represents all the data required to pass the TESTS */

/* this is a QA/QC file that is generated for the model     */
/* author to verify that the results are repeatable         */
/* write the csv outputs that will be used to verify the    */
/* results of the variant against the author's outputs      */
/* this file represents a basic QA/QC test for the model    */
/* so that we can use it for repeatable research            */

/* for now, all the units are in imperial units             */


/* this functions writes/updates a comma separated values (CSV) file    */
/* that can be read into any common spreadsheet and given the           */
/* model equations (described in model_[variant].c), the coefficients   */
/* found in coeffs_[variant].c, the plot data (plots.variant.txt), and  */
/* the plants_[variant].txt, species_[variant].txt, these values should */
/* match those that can be computed using a hand calculator, slide rule */
/* or abbucus by the model developer                                    */
/* this should be the file that is required for payment.                */ 

/* species array must be in numerical order                             */
void write_cips_test_file(
    FILE                    *fp,
    long                    cycle,

    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,

    unsigned long           n_points,
    struct PLOT_RECORD      *plots_ptr,

    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    
    unsigned long           n_years_after_planting )

{

    unsigned long           i;
    unsigned long           j;
    unsigned long           k;
    struct  PLOT_RECORD     *plot_ptr;
    struct  PLANT_RECORD    *plant_ptr;

    /* debugging/compliance labels */

//#ifdef _DEBUG
     /* these labels are used for the compliance/debugging file gen     */
    char            INT_SPARE_LABEL[40][32];
    char            DBL_SPARE_LABEL[40][32];
//#endif


    for( k = 0; k < 40; k++ )
    {
        /* strset is not ANSI so don't use it. */
        //strset( INT_SPARE_LABEL[k], '\0' );
        //strset( DBL_SPARE_LABEL[k], '\0' );

        memset( INT_SPARE_LABEL[k], '\0', 32 );
        memset( DBL_SPARE_LABEL[k], '\0', 32 );

    }


    /* assign the debugging labels */
    strcpy( DBL_SPARE_LABEL[0], "rel_height" );
    strcpy( DBL_SPARE_LABEL[1], "rh_mod_hg" );
    strcpy( DBL_SPARE_LABEL[2], "cv_mod_hg" );
    strcpy( DBL_SPARE_LABEL[3], "delta_h_pot" );
    strcpy( DBL_SPARE_LABEL[4], "psi" );
    strcpy( DBL_SPARE_LABEL[5], "base_model" );
    strcpy( DBL_SPARE_LABEL[6], "cv_mod_dg" );
    strcpy( DBL_SPARE_LABEL[7], "rh_mod_dg" );
    strcpy( DBL_SPARE_LABEL[8], "prob_of_mortality" );
    strcpy( DBL_SPARE_LABEL[9], "cw_growth" );


    /* assign the integer spares and debugging labels */
    //plant_ptr->INT_SPARE[0]         =   after_first_season;

    /* assign the debugging labels */
    strcpy( INT_SPARE_LABEL[0], "after_first_season" );


    /* sort by plot, plant */
    /* sort the tree list just to be safe   */
    qsort(   plants_ptr, 
       n_plants, 
       sizeof( struct PLANT_RECORD ), 
       compare_plants_by_plot_plant ); 

    /* order the plants by plot and plant */
            /* print the cycle number */
            fprintf( fp, "cycle," );

            /* print the headers for the plot data */
            fprintf( fp, 
                    "plot,"                      /* plot id                         */
                    "site_30,"				    /* site index, base age 30 years */
                    "error,"			            /* error flag                      */    
                    "shrub_pct_cover,"           /* percent cover for shrubs        */
                    "shrub_mean_height,"         /*  temp variables                  */
                    "shrub_expf,"                /*  shrub expansion factor          */
                    "basal_area,"                /*  total stand ba at 4.5 feet      */
                    "d6_area,"                   /*  basal area at 6inch             */
                    "expf,"                      /*  total trees per acre            */
                    "bh_expf,"                   /*  expf in trees hw+con over 4.5 ft*/
                    "d12ba_c," );				    /* like the d6_area which is the	*/

            fprintf( fp, 
                    "plot,"            
                    "plant,"           
                    "sp_code,"         
                    "d6,"              
                    "d6_area," 
                    "d12," 
                    
                    "dbh,"             
                    "basal_area,"      
                    "tht,"             
                    "cr,"			    
                    "n_stems,"         
                    "expf,"            
                    "crown_width,"     
                    "crown_area,"     
                    "user_code," );       
            
            /* put the variant specific variables in here           */
            /* these will come from the specific process            */
            /* between the model developer and source maintainer    */
            for( k = 0; k < 40; k++ )
            {

                if( strcmp( "", INT_SPARE_LABEL[k] ) != 0 )
                {
                    fprintf( fp, "%s,", INT_SPARE_LABEL[k] );
                }

            }


            for( k = 0; k < 40; k++ )
            {
                if( strcmp( "", DBL_SPARE_LABEL[k] ) != 0 )
                {
                    fprintf( fp, "%s,", DBL_SPARE_LABEL[k] );
                }

            }



            /* these are the variables that must be validated against */
            fprintf( fp, 
                    "d6_growth,"       
                    "dbh_growth,"      
                    "tht_growth,"      
                    "cr_growth,"       
                    "cw_growth,"       
                    "expf_change" );

            fprintf( fp, "\n" );



    /* for each plot, write out the plot level variables,   */
    /* this value will be repeated once for each tree       */
    /* for each plot, project it forward one year */
    plot_ptr = &plots_ptr[0];
    for( i = 0; i < n_points; i++, plot_ptr++ )
    {
        /* loop over the plants and write out a record */
        plant_ptr = &plants_ptr[0];
        for( j = 0; j < n_plants; j++, plant_ptr++ )
        {
            /* if the plant is on the plot, then    */
            /* write out the QA/QC record to the    */
            /* file                                 */
            if( plant_ptr->plot == plot_ptr->plot )
            {

                /* print out which cycle you're on */
                fprintf( fp, "%ld,", cycle );

                /* print the plot level data to the QA/QC test file */
                fprintf( fp, 
                        "%ld,%lf,%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,",
                        plot_ptr->plot,                      /* plot id                         */
                        plot_ptr->site_30,				    /* site index, base age 30 years */
                        plot_ptr->error,			            /* error flag                      */    
                        plot_ptr->shrub_pct_cover,           /* percent cover for shrubs        */
                        plot_ptr->shrub_mean_height,         /*  temp variables                  */
                        plot_ptr->shrub_expf,                /*  shrub expansion factor          */
                        plot_ptr->basal_area,                /*  total stand ba at 4.5 feet      */
                        plot_ptr->d6_area,                   /*  basal area at 6inch             */
                        plot_ptr->expf,                      /*  total trees per acre            */
                        plot_ptr->bh_expf,                   /*  expf in trees hw+con over 4.5 ft*/
                        plot_ptr->d12ba_c );				    /* like the d6_area which is the	*/
//                        plot_ptr->d12ba_c * FTAC2M2HA );				    /* like the d6_area which is the	*/
                                                                /* basal area for "conifers"		*/
                                                                /* (aka DF) at 30 cm (15 inches)	*/
                                                                /* above the ground					*/

                fprintf( fp, 
                        "%ld,%ld,%s,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%ld,%lf,%lf,%lf,%ld,",
                        plant_ptr->plot,            
                        plant_ptr->plant,     
                        
                        species_ptr[plant_ptr->sp_idx].sp_code,
                        
                        plant_ptr->d6,              
                        plant_ptr->d6_area,
                        plant_ptr->d12,
                        plant_ptr->dbh,             
                        plant_ptr->basal_area,      
                        plant_ptr->tht,             
                        plant_ptr->cr,			    
                        plant_ptr->n_stems,         
                        plant_ptr->expf,            
                        plant_ptr->crown_width,     
                        plant_ptr->crown_area,      
                        plant_ptr->user_code );   


                /* print the variant specific (intermediate) values here    */
                /* put the variant specific variables in here               */
                /* these will come from the specific process                */
                /* between the model developer and source maintainer        */
                //for( k = 0; k < N_PLANT_INT_SPARES_LEFT; k++ )
                for( k = 0; k < 40; k++ )
                {
                    if( strcmp( "", INT_SPARE_LABEL[k] ) != 0 )
                    {
                        fprintf( fp, "%ld,", plant_ptr->INT_SPARE[k] );
                    }

                }


                
                for( k = 0; k < 40; k++ )
                {
                    if( strcmp( "", DBL_SPARE_LABEL[k] ) != 0 )
                    {
                        fprintf( fp, "%.10lf,", plant_ptr->DBL_SPARE[k] );
                    }
                }

                /* at the end of the projection, plot the increment variables   */
                /* these are what the modeler needed to generate                */
                fprintf( fp, 
                        "%lf,%lf,%lf,%lf,%lf,%lf,",
                        plant_ptr->d6_growth,       
                        plant_ptr->dbh_growth,      
                        plant_ptr->tht_growth,      
                        plant_ptr->cr_growth,       
                        plant_ptr->cw_growth,       
                        plant_ptr->expf_change );

                fprintf( fp, "\n" );


            } /* end of if-then statement for plants on plot */

        } /* end of loop for plants */
   
    } /* end of loop for the plots */



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




/* include this modification today */


/*
Equation for CIPS-CONIFERS; for calculating D6 and D12

Diameter at any height = 0.4361*(Vegcov– 0.0504) * (H0.8046)*(X (-11.8347*(h/HT)+(5.7161*(h/H)**2) 
    - (0.2206*ln((h/H)+0.0001)) + (11.3148*(h/H)**0.5) - (1.9433*exp(h/H))+ (6.3261*I*(dbh/H)))

where	dbh		=	dbh (mm)
	H		=	Height (cm)
	h		=	Arbitrary height (cm) on stem with 0=h=H
	Vegcov	=	Competing vegetation cover (%)
	X		= 	2*(1-(h/H)^0.5)
	I		=	1 if H>137 cm; 0 otherwise

*/


/********************************************************************************/
/*                  calc_dhi                                                    */
/********************************************************************************/
/*  Description :                                                               */
/* this function computes the basal diameter for the tree record, given the     */
/* total height and the proportion of the area covered by competing vegetation  */
/* on the plot                                                                  */
/*  Author      :   Jeff D. Hamann and Douglas A. Maguire                       */
/*  Date        :   June 20, 2011                                               */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  pct_veg_cover               -   total height of the subject tree            */
/*  *pred_dbh                    -   predicted breast height diameter           */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : Not in Annual Report Document                                     */
/*  Source  : Center for Intensive Planted-Forest Silviculture,                 */
/*              2010 Annual Report                                              */
/*              http://www.fsl.orst.edu/cips/                                   */
/*  Coeffs  : d6 (repurposed, actually same equation)                           */
/********************************************************************************/
void cips_calc_dob_hi( 
    unsigned long           *return_code,
    struct PLOT_RECORD      *plot_ptr,
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  hi,
    double                  *di )
{

    double  part_a  = 0.0;
    //double  part_b  = 0.0;  //unused removed jan 2014 mwr;
    double  X       = 0.0;
    double  hi_tht  = 0.0;

    double  b1      = 0.0;
    double  b2      = 0.0;
    double  b3      = 0.0;
    double  b4      = 0.0;
    double  b5      = 0.0;
    double  b6      = 0.0;
    double  b7      = 0.0;
    double  b8      = 0.0;
    double  b9      = 0.0;
    
    unsigned long I =   0;

    /* initialize the return diameter to zero */
    *di = 0.0;

    /* perform error check for the correct number of coeffs */
    if( coeffs_ptr == NULL )
    {
        *return_code = INVALID_COEFF;
        return;
    }

    if( plant_ptr->tht * FT2CM <= 15.0 )
    {
        *return_code = INVALID_INPUT_VAL;
        return;
    }

    //b1 = 0.4361;
    //b2 = -0.0504;
    //b3 = 0.8046;
    //b4 = -11.8347;
    //b5 = 5.7161;
    //b6 = 0.2206;
    //b7 = 11.3148;
    //b8 = 1.9433;
    //b9 = 6.3261;

    b1  = coeffs_ptr->dob_hi[0];  /*  0.4361;     */
    b2  = coeffs_ptr->dob_hi[1];  /* -0.0504;     */
    b3  = coeffs_ptr->dob_hi[2];  /* 0.8046;      */
    b4  = coeffs_ptr->dob_hi[3];  /* -11.8347;    */
    b5  = coeffs_ptr->dob_hi[4];  /* 5.7161;      */
    b6  = coeffs_ptr->dob_hi[5];  /* 0.2206;      */
    b7  = coeffs_ptr->dob_hi[6];  /* 11.3148;     */
    b8  = coeffs_ptr->dob_hi[7];  /* 1.9433;      */
    b9  = coeffs_ptr->dob_hi[8];  /* 6.3261;      */

    /* to prevent raising zero to a power */
    /* add a little tiny value to the pct_veg_cover if it's zero */
    if( plot_ptr->shrub_pct_cover <= 0.000001f )
    {
        plot_ptr->shrub_pct_cover = 0.000001f;
    }

    /* now equation converted from the equation emailed from Doug Mainwaring */
//    where	dbh		=	dbh (mm)
//	H		=	Height (cm)
//	h		=	Arbitrary height (cm) on stem with 0=h=H
//	Vegcov	=	Competing vegetation cover (%)

//	X		= 	2*(1-(h/H)^0.5)
    X		= 	2.0 * ( 1.0 - pow( ( hi * FT2CM ) / ( plant_ptr->tht * FT2CM ), 0.5 ) );

//	I		=	1 if H>137 cm; 0 otherwise
    I = 0;
    if( plant_ptr->tht * FT2CM > 137.0f )
    {
        I = 1;
    }

//Diameter at any height = 0.4361*(Vegcov– 0.0504) * (H0.8046)*(X (-11.8347*(h/HT)+(5.7161*(h/H)**2) 
//    - (0.2206*ln((h/H)+0.0001)) + (11.3148*(h/H)**0.5) - (1.9433*exp(h/H))+ (6.3261*I*(dbh/H)))

/*
    Diameter at any height = 0.4361 *   ( Vegcov– 0.0504 ) * 
                                    ( H0.8046 ) *
                                    ( X (-11.8347 * (h/HT)+
                                    ( 5.7161 * (h/H)**2) -
                                    ( 0.2206 * ln((h/H)+0.0001)) + 
                                    ( 11.3148 * (h/H)**0.5) - 
                                    ( 1.9433 * exp(h/H)) + (6.3261*I*(dbh/H) ) )
*/


    hi_tht = ( hi * FT2CM ) / ( plant_ptr->tht * FT2CM );

/*
    part_a = -11.8347 * hi_tht +
               5.7161 * pow( hi_tht, 2.0 ) -
               0.2206 * log( hi_tht + 0.0001 ) + 
              11.3148 * pow( hi_tht, 0.5 ) - 
               1.9433 * exp( hi_tht ) + 
               6.3261 * I * ( plant_ptr->dbh / plant_ptr->tht * FT2CM );
*/

    /* todo: have to add some code here to make sure this doesn't return    */
    /* bad numbers.                                                         */
    part_a = b4 * hi_tht +
             b5 * pow( hi_tht, 2.0 ) -
             b6 * log( hi_tht + 0.0001 ) + 
             b7 * pow( hi_tht, 0.5 ) - 
             b8 * exp( hi_tht ) + 
             b9 * I * ( ( plant_ptr->dbh * IN2MM ) / ( plant_ptr->tht * FT2CM ) );

    /* entered from contract */
    *di = b1 * pow( plot_ptr->shrub_pct_cover, b2 ) * 
                pow( ( plant_ptr->tht * FT2CM ), b3 ) * 
                pow( X, part_a );

    /* the function returns in inches */
    *di *= MM2IN;
    *return_code = CONIFERS_SUCCESS;

    if( *di < 0.0f )
    {
        *di   = 0.0f;
        *return_code    = CONIFERS_ERROR;
    }


}


