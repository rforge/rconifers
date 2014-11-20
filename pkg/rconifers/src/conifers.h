/****************************************************************************/
/*                                                                          */
/*  conifers.h                                                              */
/*  conifers header file                                                    */
/*                                                                          */
/****************************************************************************/
                                                                                
/* 	$Id: conifers.h 880 2012-04-10 22:35:46Z hamannj $	 */

#ifndef __CONIFERS_H__
#define __CONIFERS_H__

/* this only applies if you're on windows */
#ifndef WIN32
#define __stdcall
#endif

#ifdef __cplusplus  
extern "C" {
#endif

/****************************************************************************/
/* Definitions                                                              */
/****************************************************************************/
#define SP_LENGTH           6           /*  max length of sp_codes          */   
#define FVS_SP_LENGTH       3           /*  max length of sp_codes for fvs  */   
#define NAME_LENGTH         40          /*  max length for long names       */   
#define MAX_COEFFS          16          /*  current max coeffs for dynamic  */
#define MIN_HEIGHT          0.5         /* minimum total height             */
#define	AIT_BIN_RES         0.10        /* single array for bal/cat         */
#define	AIT_SIZE            750         /* single array for bal/cat         */

#ifndef FALSE                           /* if FALSE is not already defined  */
#define FALSE               0
#endif

#ifndef TRUE                            /* if TRUE is not already defined   */
#define TRUE                1
#endif

#define FC_I        5.454153912482e-3   /*  forestry constant imperial      */
#define FC_M	    7.85398163398e-5    /*  forestry constnt metric         */ 
#define MY_PI       3.14159265359       /*                                  */    
#define ONE_OVER_PI 0.318309886184      /*                                  */   
#define DEG2RAD     1.74532925199e-2    /*  degrees to radians              */   

#define FT2CM     30.48000		        /*  feet to centimeters             */
#define CM2FT     1.0/FT2CM		        /*  centimeters to feet		        */

#define FT2M      0.304800              /*  feet to meters                  */
#define GRM2LB    0.0022046226          /*  grams to pounds                 */
#define KG2LB     2.2046226             /*  kilograms to pounds             */
#define LB2TON    0.00050               /*  pounds to tons                  */
#define IN2CM     2.54000				/*  inches to centimeters			*/
#define IN2MM     25.4000				/*  inches to millimeters           */
#define MM2IN	  1.0/IN2MM			    /*  millimeters to inches           */
#define FTAC2M2HA 0.229568411			/* ft^2/acre to m^2/ha	  		    */
#define M2HA2FTAC 1.0/FTAC2MHA		    /* ft^2/acre to m^2/ha		        */


#define REINEKE_B1  1.605               /*  Reineke's constant              */   
#define REINEKE_B1A 1.77                /*  Modified by Oliver et al.       */   
#define CRITICAL_RD 0.600               /*  threshold imminent sdi mort     */   
#define CCF_CONST_I 0.00180302608677    /*  PI / 4.0 / 43560.0 * 100.0      */   
#define SQ_FT_PER_ACRE  43560.0         /*  Square feet per acre            */   

/* todo: these varlues need to move into the coeffs */
//#define CURRENT_COEFFS_VER  4.14        /*  current version num for coeffs  */
//#define MODEL_VERSION       4.14        /* model (equations) version        */

/* error constants          */
#define CONIFERS_SUCCESS        0       
#define CONIFERS_ERROR          1
#define INVALID_COEFF           2
#define INVALID_EQUATION        3       
#define INVALID_FSP             4
#define INVALID_OPTION          5
#define INVALID_INPUT_VAL       6
#define INVALID_SP_CODE         7
#define INVALID_FILE_NAME       8
#define CONIFERS_MORT_ERROR     9
#define THINNING_ERROR          10
#define BELOW_THIN_TARGET       11  /*  used in thin.c plot is below target     */
#define WRONG_COEFF_VERSION     12  /*  MOD043a coeffs version check            */
#define FILL_SAMPLE_FAILED      13  /*  MOD050 failed to fill in awi            */
#define SYS_ARCHIVE_FAILED      14  /*  MOD051 systum archive import failed     */
#define SAMPLE_RECORD_FAILED    15  /*  MOD051 failed to create sample record   */
#define CORRUPTED_ARCHIVE       16  /*  MOD051 bad systum-1 archive file        */
#define INVALID_PLANT_COUNT     17  /*  MOD068 */
#define INVALID_PLOT_COUNT      18  /*  MOD068 */
#define INVALID_BIN_COEFFS_FILE 19  
#define INVALID_VARIANT		    20  
#define INVALID_PLANT_TYPE      21
#define FAILED_MEMORY_ALLOC     22
#define FAILED_PROJECT_PLANT    23
#define INVALID_SPECIES         24
#define INVALID_DBH             25
#define INVALID_HEIGHT          26
#define FILL_VALUES_ERROR       27


/* plant record error constants (hex values)            */
#define E_OKDOKEY                 0x00000000  /*    0 */
#define E_INVALID_SPECIES         0x00000001  /*    1 */
#define E_INVALID_DBH             0x00000002  /*    2 */
#define E_INVALID_HEIGHT          0x00000004  /*    4 */
#define E_INVALID_EXPF            0x00000008  /*    8 */
#define E_INVALID_CR              0x00000010  /*   16 */
#define E_INVALID_PCT_COVER       0x00000020  /*   32 */
#define E_INVALID_CW              0x00000040  /*   64 */
#define E_INVALID_MCW             0x00000080  /*  128 */
#define E_INVALID_HG              0x00000100  /*  256 */
#define E_INVALID_D6G             0x00000200  /*  512 */
#define E_INVALID_DBHG            0x00000400  /* 1024 */
#define E_INVALID_CRG             0x00000800  /* 2048 */
#define E_INVALID_D6              0x00001000  /* 4096 */

/* thinning constants */
#define DO_SDI_MORT                 0   /*  thinning type 0 is sdi mort     */
#define DO_EXPF_SP_THIN             1   /*  thin a single sp on expf        */
#define DO_EXPF_THIN                2   /*  thin a single sp on expf        */
#define DO_EXPF_THIN_FROM_BELOW     3   /*  thin from below based on dbh    */
#define DO_EXPF_SP_THIN_FROM_BELOW  4   /*  thin from below based on dbh-sp */

/* activity switches for sorting the summaries */
#define NO_ACTIVITY                 0
#define THIN_ACTIVITY               1
#define SHRUB_CONTROL_ACTIVITY      2

/* added for conifers 3.0   */
#define CONIFER					0
#define HARDWOOD                1
#define SHRUB                   2
#define FORB                    3
#define NON_STOCKED             4  
#define PLANT_TYPES             5

/* variants added for conifers 4.0 */
/* todo: step #1 - add new variant #define here */
#define CONIFERS_SWO            0
#define CONIFERS_SMC            1
#define CONIFERS_SWOHYBRID      2
#define CONIFERS_CIPS           3
// #define N_VARIANTS              (CONIFERS_CIPS+1)

/*
variant_id,
coeffs_version, 
model_version,
citation,
*/


/* We might need a structure to hold the contents   */
/* of the simulator variant                         */
/* possible data required might include             */
/* id, name, metric/imperial, n_species             */

/* todo: step #2 - create new coeffs_[VARIANT].c file, enter, and verify the coefficients */
/* todo: step #3 - create new model_[VARIANT].c file, code, and verify the functions  */




/****************************************************************************/
/* Structure Definitions                                                    */
/****************************************************************************/

/* The members for this structure are contain the coeffs that       */
/* can be read in from a text or binary file using the functions    */
/* found in file_io.c. The nomenclature for the coeffs names        */
/* represent the predicted variable followed by the explanitory     */
/* variables. For example dbh_d6_height, represents the coeffs      */
/* for the equation that predict dbh from d6 and total height       */
/* For now, there are 20 possible equations that can be saved       */
/* for each functional species coeffs set                           */
   
   struct COEFFS_RECORD 
   {
	 unsigned long   idx;        /* this is the functional species code -- group in 2.0 */
	 unsigned long   type;       /* CONIFER, HARDWOOD, SHRUB, or FORB */

	 /* this is required for the user interface */
	 char           group[SP_LENGTH];               /*  name abbreviated from the usfs  */

	 double         d6_growth[MAX_COEFFS];          /*  D1: growth coeffs               */
	 double         ht_growth[MAX_COEFFS];          /*  D3: growth coeffs               */ 
	 double         cr_growth[MAX_COEFFS];          /*  D4: hcb change parameters       */
	 double         crown_ratio[MAX_COEFFS];      /*  S3: crown ratio coeffs          */
	 double         crown_width[MAX_COEFFS];      /*  S1: crown width coeffs, old CA  */
	 double         max_crown_width[MAX_COEFFS];  /*  S2: crown width coeffs, old CA  */
	 double         d6_ht_dbh[MAX_COEFFS];        /*  S5, S7                          */
	 double         d6_ht[MAX_COEFFS];            /*  S4, S9                          */
	 double         dbh_ht[MAX_COEFFS];           /*  S6, S10                         */
	 double         n_stems_from_ht[MAX_COEFFS];  /*  S11                             */
	 double         dbh_growth[MAX_COEFFS];       /*  XX: dbh growth coefficients     */

	 double         cfvolume4[MAX_COEFFS];        /*  S12 volume coefficients         */
	 double         biomass[MAX_COEFFS];          /*  S13 biomass coefficients        */
	 double         cw_growth[MAX_COEFFS];          /*  d6 cw_growth coefficients        */
   
    /* coefficients added for the CIPS variant */
    /* todo: you need to added these the old models too */
    double	        dbh_ht_veg_cov[MAX_COEFFS]; /* dbh-ht prediction            */
    double	        dob_hi[MAX_COEFFS];         /* dob at height hi             */	
    double          mortality[MAX_COEFFS];      /* species specific mortality   */

	 /* add new coefficients here */
     /* make sure you include the copy in the variant specific files */
   
   };


/* this structure maintains the relationship between the                    */
/* actual species entered by the user and the functional species that       */
/* drives the model. See functions get_species_entry and get_coeffs_entry   */
   struct SPECIES_RECORD
   {
      unsigned long  idx;
      unsigned long  fsp_idx;
      
      char           sp_code[SP_LENGTH];             /*   usfs species code              */
      char           common_name[NAME_LENGTH];       /*  common plant name               */
      long           organon_sp_code;                /*  used to translate to organon    */
      long           cactos_sp_code;                 /*  used to translate to cactos     */
      char           fvs_sp_code[FVS_SP_LENGTH];     /*  used to translate to fvs        */
      
      /*  mortality coefficients          */
      double		endemic_mortality;			     /*  times 100     */
      double		max_sdi;        
      double		browse_damage;
      double		mechanical_damage;

      /* genetic coeffs here */
      double		genetic_worth_h;
      double		genetic_worth_d;

      /* climate model coeffs	*/
      double		min_temp;
      double		max_temp;
      double		opt_temp;
      
   };


   struct PLANT_RECORD 
   {

	 unsigned long  plot;            /*  plot id                         */
	 unsigned long  plant;           /*  tree id                         */
	 unsigned long  sp_idx;          /*  index into the species array    */
	 double	        d6;              /*  diameter at 6"                  */
	 double         d6_area;         /*  area (ft^2) of d6 cross section */
	 double	        dbh;             /*  diameter at breast height       */
	 double         basal_area;      /*  dbh basal area                  */
	 double	        tht;             /*  total height                    */
	 double	        cr;			     /*  crown ratio (a proportion value)*/
	 unsigned long	n_stems;         /*  number of stems                 */
	 double	        expf;            /*  expansion factor                */
	 double         pct_cover;       /*  percent cover                   */
	 double         crown_width;     /*  root of crown_area*4/pi         */
	 double         crown_area;      /*  crown area @ max width          */
	 unsigned long  user_code;       /*  generic user code               */
	 double         d6_growth;       /*  change in d6                    */
	 double         dbh_growth;      /*  change in dbh                   */
	 double         tht_growth;      /*  1-year total height growth      */
	 double         cr_growth;       /*  change in crown ratio           */
	 double         cw_growth;       /*  change in crown width           */
	 double         expf_change;     /*  number of trees killed          */
	 unsigned long  errors;          /*  simple error flag for the tree  */
	 double         max_crown_width; /* maximum cw from Paine&Hann       */

     /* added for the CIPS model */
	 double         d12;			 /* diameter at 30 cm				*/
	 double         d12_growth;		 /* change in diameter at 30 cm		*/
	 double         d12_area;		 /* area (ft^2) of d12 cross section */

	 /* spares. these are reserved for debugging, new variables, etc.    */
	 long           INT_SPARE[30];   /*  generic spares, reserved        */
	 double         DBL_SPARE[30];	 /*  for debugging variants          */

   };

   struct PLOT_RECORD
   {

	 unsigned  long plot;                   /*  plot id                         */
	 double         latitude;               /*  plot latitude, in dd MOD010     */
	 double         longitude;              /*  plot longitude, in dd MOD010    */
	 double         elevation;              /*  plot elevation above msl, ft    */
	 double         slope;                  /*  slope in percent?               */
	 double         aspect;                 /*  aspect ( 0 - 360 degrees )      */
	 double         water_capacity;         /*  water holding capactiy          */
	 double         mean_annual_precip;     /*  annual precip                   */
	 double		    site_30;				/* site index, base age 30 years */
	 unsigned long	error;			/*  error flag                      */    
	 double         shrub_pct_cover;        /*  percent cover for shrubs        */
	 double         shrub_mean_height;      /*  temp variables                  */
	 double         shrub_expf;             /*  shrub expansion factor          */
	 double         basal_area;             /*  total stand ba at 4.5 feet      */
	 double         d6_area;                /*  basal area at 6inch             */
	 double         expf;                   /*  total trees per acre            */
	 double         hann_wang_x0;           /*  hann & wang x0 value mort init  */
	 double         bh_expf;                /*  expf in trees hw+con over 4.5 ft*/
	 double         qmd;                    /*  quadratic mean diameter         */
	 double         sdi;                    /*  stand density index             */
	 double         sdimax;                 /*  maximum stand dens index        */
	 double         crown_area;             /*  crown area sq ft                */
	 double         ccf;                    /*  crown competition factor        */
	 double         ba_c;                   /*  b.h. basal area in conifers     */
	 double         ba_h;                   /*  b.h. basal area in hardwoods    */
	 double         d6ba_c;                 /*  basal ba in conifers            */
	 double         d6ba_h;                 /*  basal ba in hardwoods           */
	 double         d6ba_s;                 /*  basal area in shrubs            */                                    
	 double         ca_c;                   /*  crown area in conifers          */
	 double         ca_h;                   /*  crown area in hardwoods         */
	 double         ca_s;                   /*  crown area in shrubs            */
	 double         bait[PLANT_TYPES][AIT_SIZE];  /* basal area in taller       */
	 double         cait[PLANT_TYPES][AIT_SIZE];  /* crown area in taller       */
	 

	 /* variables added for the CONIFERS_CIPS model */
	 double			d12ba_c;				/* like the d6_area which is the	*/
											/* basal area for "conifers"		*/
											/* (aka DF) at 30 cm (15 inches)	*/
											/* above the ground					*/
	 double         bal[PLANT_TYPES][AIT_SIZE];  /* basal area in larger */
                                                 /* diameter trees */

	 /* variables added for the CONIFERS_SWOHYBRID model */
     double			growing_season_precip;  /* this is not mean annual precip */
	 double			solar_radiation[12];	/*	solar radiation, MegaJoules/m^2	*/
	 double			mean_monthly_temp[12];	/*	mean monthly temp in C			*/

	 /* spares. these are reserved for debugging, new variables, etc.    */
	 long           INT_SPARE[30];   /*  generic spares, reserved        */
	 double         DBL_SPARE[30];	 /*  for debugging variants          */

   };

/* This structure serves as a general purpose structure to hold         */
/* calculated statistics used in reporting, growth, and diagnostics     */
/* this structure could be embedded into the PLOT_RECORD becuase        */
/* the variables that are required in the PLOT_RECORD are also present  */
/* here.                                                                */
   struct SUMMARY_RECORD
   {
	 long           time;                   /*  time or age for the summaries   */
	 long           activity;               /*  activity switch for sorting     */
	 
	 unsigned long  code;                   /*  sp_code, fsp, type, anything    */
	 
	 double         min_dbh;                /*  minumim dbh                     */
	 double         qmd;                    /*  quad mean diameter              */
	 double         max_dbh;                /*  max dbh in the list             */
	 double         basal_area;             /*  basal area for the code         */
	 double         min_height;             /*  min total height                */
	 double         mean_height;            /*  average total height            */
	 double         max_height;             /*  max total height                */
	 double         cr;                     /*                                  */
	 double         expf;                   /*  total stand stems per acre      */
	 double         bh_expf;                /*  expf for hw+con over 4.5 feet   */
	 double         tree_expf;              /*  expf for trees only             */
	 double         crown_area;             /*                                  */
	 double         ccf;                    /*  crown comp factor               */
	 double         pct_cover;              /*                                  */
	 double         sdi;                    /*  total stand sdi                 */
	 double         sdimax;                 /*  stand maxsdi                    */
	 double         rel_density;            /*                                  */
	 double         curtis_rd;              /*  curtis relative density         */
	 double         cfvolume4;              /*  cf volume to a 4 " top          */
	 double         biomass;                /*  biomass                         */

	 double         min_hd6_ratio;          /*  min ht/d6 ratio                 */
	 double         mean_hd6_ratio;         /*  average ht/d6 ratio             */
	 double         max_hd6_ratio;          /*  max ht/d6 ratio                 */
	 double         height_40;              /*                                  */
	 double         con_tpa;                /*  conifers trees per acre         */

     /* spares   */
	 long           INT_SPARE[30];          /*  generic spares, reserved        */
	 double         DBL_SPARE[30];          /*  generic spares, reserved        */


   };




/****************************************************************************/
/* function declarations                                                    */
/****************************************************************************/


/****************************************************************************/
/* functions in coeffs_mgt.c                                                */
/****************************************************************************/
struct COEFFS_RECORD *con_init_coeffs( 
	unsigned long				version,
	unsigned long				*n_coeffs,
    double                      *coeffs_version,
    double                      *model_version );

   struct COEFFS_RECORD *con_swo_init_coeffs( 
      unsigned long         *n_coeffs, 
    double                      *coeffs_version,
    double                      *model_version );

   struct COEFFS_RECORD *con_smc_init_coeffs( 
      unsigned long         *n_coeffs, 
    double                      *coeffs_version,
    double                      *model_version );

   struct COEFFS_RECORD *con_cips_init_coeffs( 
      unsigned long         *n_coeffs, 
    double                      *coeffs_version,
    double                      *model_version );

struct COEFFS_RECORD *con_swo_hybrid_init_coeffs( 
   unsigned long               *n_coeffs, 
    double                      *coeffs_version,
    double                      *model_version );

/* todo: make sure you include the function declarations for the init_coeffs */
/* functions here */


/****************************************************************************/
/* functions in conifers.c                                                  */
/****************************************************************************/
   unsigned long __stdcall is_tree( 
      struct COEFFS_RECORD *c_ptr );
   
   unsigned long __stdcall is_shrub( 
      struct COEFFS_RECORD *c_ptr );

   unsigned long __stdcall is_forb( 
      struct COEFFS_RECORD *c_ptr );

   unsigned long __stdcall is_non_stocked( 
      struct COEFFS_RECORD *c_ptr );

   void __stdcall init(	unsigned long *return_code, 
			unsigned long variant );


/****************************************************************************/
/* functions in file_io.c                                                   */
/****************************************************************************/

/* move all the internal and read functions to the top */
   struct PLANT_RECORD *combine_plant_lists( 
      unsigned long           *return_code,
      unsigned long           *n_total_plants,
      unsigned long           n_new_plants,
      struct PLANT_RECORD     *original_plants_ptr,
      struct PLANT_RECORD     *new_plants_ptr );
   
   void free_sample_data( 
      unsigned long             *n_points,
      struct PLOT_RECORD        **plots_ptr,
      unsigned long             *n_plants,
      struct PLANT_RECORD       **plants_ptr );
   

   void read_systum1_file( 
      unsigned long           *return_code,
      const char              *filename, 
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           *n_points,
      struct PLOT_RECORD      **plots_ptr,
      unsigned long           *n_plants,
      struct PLANT_RECORD     **plants_ptr,
      unsigned long			  *age );
	

   void __stdcall dump_plots_to_file( 
      unsigned long           *return_code,
      const char              *filename,
      unsigned    long        n_points,
      struct PLOT_RECORD      *plots_ptr );

   void __stdcall dump_plants_to_file( 
      unsigned long			  *return_code,
    const char          *filename,
//    FILE                *fp,
      unsigned long			  n_plants,
      struct PLANT_RECORD	  *plants_ptr,
      unsigned long			  n_species,
      struct SPECIES_RECORD	  *species_ptr );

   void __stdcall dump_larger_to_file( 
      unsigned long           *return_code,
      const char              *filename,
      unsigned    long        n_points,
      struct PLOT_RECORD      *plots_ptr );

   void __stdcall write_summaries_to_file( 
      unsigned long           *return_code,
      const char              *filename,
      unsigned long           n_records,
      struct SUMMARY_RECORD   *summaries_ptr,
      unsigned long			  n_species,
      struct SPECIES_RECORD	  *species_ptr,
      unsigned long           age);

   void read_plots_from_file( 
      unsigned long			*return_code,
      const char				*filename, 
      unsigned long			*n_points,
      struct PLOT_RECORD		**plots_ptr );

   void read_plants_from_file( 
	unsigned long			*return_code,
    const char				*filename, 
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long			*n_plants,
	struct PLANT_RECORD		**plants_ptr );

   void __stdcall read_sample_from_file( 
      unsigned long			*return_code,
      const char				*filename,
      unsigned long			n_species,
      struct SPECIES_RECORD	*species_ptr,
      unsigned long           *n_points,
      struct PLOT_RECORD      **plots_ptr,
      unsigned long           *n_plants,
      struct PLANT_RECORD     **plants_ptr,
      unsigned long			*age,
      unsigned long         *variant);

   void __stdcall write_plots_to_text_file( 
      unsigned long       *return_code,
      const char          *filename, 
      unsigned long       n_records,
      struct PLOT_RECORD  *plots_ptr );

   void __stdcall write_sample_to_file( 
      unsigned long           *return_code,
      const char              *filename, 
      unsigned long			  age,
      unsigned long           variant,
      unsigned long           n_points,
      struct PLOT_RECORD      *plots_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr );

   void __stdcall write_organon_file( 
      unsigned long           *return_code,
      const char              *filename, 
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      struct COEFFS_RECORD    *coeffs_ptr);

   void __stdcall write_cactos_file( 
      unsigned long           *return_code,
      const char              *filename, 
      const char              *sample_id, 
      unsigned long           n_points,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_site_indecies,
      double                  *site_indecies_ptr,     
      unsigned long           n_ages,
      double                  *ages_ptr,              
      unsigned long           n_coeffs,               
      struct COEFFS_RECORD    *coeffs_ptr );          

   void __stdcall write_cactos_ingrowth_file(
      unsigned long           *return_code,
      const char              *filename, 
      unsigned long           n_points,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,               /* MOD056 */
      struct COEFFS_RECORD    *coeffs_ptr );          /* MOD056 */

   void __stdcall write_fvs_file( 
      unsigned long           *return_code,
      const char              *filename, 
      unsigned long           variant,
      /* const char          *stand_id,       */
      /* double              stand_age,       */
      /* double              current_year,    */
      unsigned long           n_points,
      struct PLOT_RECORD      *plots_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr    );

   void __stdcall write_systum1_file( 
      unsigned long			*return_code,
      const char				*filename, 
      unsigned long			n_trees,
      struct PLANT_RECORD		*tree_ptr,
      unsigned long			n_species,
      struct SPECIES_RECORD	*species_ptr );

//struct SAMPLE_RECORD *read_systum1_archive(
//    unsigned long           *return_code,
//    const char              *filename, 
//    unsigned long           n_species,
//    struct SPECIES_RECORD   *species_ptr );

   void read_systum1_archive(
      unsigned long           *return_code,
      const char              *filename, 
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,	
      unsigned long			*age,
      unsigned long           *n_points,
      struct PLOT_RECORD      **plots_ptr,
      unsigned long           *n_plants,
      struct PLANT_RECORD     **plants_ptr );

   void __stdcall write_organon_inp( 
      unsigned long           *return_code,
      const char              *filename, 
      const char              *title,             /* MOD076 */
      unsigned long           sim_version,
      unsigned long           even_age,
      unsigned long           bh_age,
      unsigned long           stand_age,
      double                  df_site,
      double                  pp_site,
      unsigned long           n_plots,
      struct PLOT_RECORD      *plots_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr    );


    void print_errors_and_warnings(
        unsigned long           n_plants,
        struct PLANT_RECORD     *plants_ptr );


/****************************************************************************/
/* functions in grow.c                                                      */
/****************************************************************************/

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
   unsigned long            *n_years_projected );

void get_taller_attribs( 
    double                  height,
    struct PLOT_RECORD      *plot_ptr,
    double                  *bait, 
    double                  *cait );

void get_age_cut(
    unsigned long           *return_code,
    double                  h40,                 
    double                  si_30,
    unsigned long           *genetics_age_cut);

/****************************************************************************/
/* functions in swo_model.c                                                 */
/****************************************************************************/
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
	double                  baf );

void swo_project_plant(  
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   struct PLANT_RECORD     *plant_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           endemic_mortality,      
   int                     hbc_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err );


/****************************************************************************/
/* functions in smc_model.c                                                 */
/****************************************************************************/
void smc_impute( 
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
	double                  baf );

void smc_project_plant(  
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   struct PLANT_RECORD     *plant_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           endemic_mortality,      
   int                     hbc_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err, 
   struct SUMMARY_RECORD   *before_sums,
   unsigned long            use_genetic_gains,
   unsigned long            genetics_age_cut);


/****************************************************************************/
/* functions in cips_model.c                                                */
/****************************************************************************/
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
	double                  baf );

void cips_project_plant(  
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   struct PLANT_RECORD     *plant_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           endemic_mortality,      
   int                     hbc_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err, 
   struct SUMMARY_RECORD   *before_sums,
   unsigned long            use_genetic_gains,
   unsigned long            genetics_age_cut,
   unsigned long			plantation_age,
   unsigned long            *n_years_projected );


/****************************************************************************/
/* functions in swo_hybrid_model.c                                          */
/****************************************************************************/
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
    double                  baf );

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
    unsigned long           use_rand_err );



/****************************************************************************/
/* functions in mortality.c                                                 */
/****************************************************************************/
   void calc_hann_wang_x0(   
      unsigned long   *return_code,
      double          sdi,
      double          bhtpa,
      double          sdimx,
      double			rd_before,
      double			rd_after,
      double          *x0 );

   void calc_sdi_mortality(
      unsigned long   *return_code,
      double          qmd,
      double          sdi,
      double          sdimax,
      double          x0,
      double          bhtpa,
      double          *mortality_proportion );

   void calc_init_x0(   
      unsigned long   *return_code,
      double          sdi,
      double          rd_current,
      double          bhtpa,
      double          sdimx,
      double          *x0 );


/****************************************************************************/
/* functions in plot.c                                                      */
/****************************************************************************/
   struct PLOT_RECORD *build_plot_array_from_plants( 
      unsigned long       *return_code,
      unsigned long       n_records, 
      struct PLANT_RECORD *plants_ptr,
      unsigned long       *n_points );


   void    copy_point_data( 
      unsigned long       *return_code,
      unsigned long       n_points,
      struct PLOT_RECORD  *src_plot_ptr,
      struct PLOT_RECORD  *dest_plots_ptr );

   void get_plant_indecies_for_plot(    
      unsigned long       *return_code,
      struct PLOT_RECORD  *plot_ptr,
      unsigned long       n_plants,
      struct PLANT_RECORD *plants_ptr,
      unsigned long       *start_idx,
      unsigned long       *end_idx,
      unsigned long       *n_plant_records_on_plot );

   struct PLOT_RECORD *get_plot( 
      long                plot,
      unsigned long       n_records,
      struct PLOT_RECORD  *plots_ptr );

   void reduce_pct_cover( 
      unsigned long           *return_code,
      double                  target_pct,
      double                  current_pct,
      unsigned long           target_sp, 
      unsigned long           n_plants,
      unsigned long           n_plots,
      struct PLANT_RECORD     *plants_ptr,
      struct PLOT_RECORD      *plots_ptr);


/****************************************************************************/
/* functions in sample.c                                                    */
/****************************************************************************/
   unsigned long calc_replication_factor(
      unsigned long           *return_code,
      unsigned long           n_points,
      unsigned long           n_plants,
      unsigned long           max_sample_size );

   struct PLOT_RECORD  *generate_duplicate_plots(
      unsigned long       *return_code, 
      unsigned long       *n_new_plots,
      unsigned long       n_dups,
      unsigned long       n_orig_plots,
      struct PLOT_RECORD  *plots_ptr );

   struct PLANT_RECORD  *generate_duplicate_plants(
      unsigned long       *return_code, 
      unsigned long       *n_new_plants,
      unsigned long       n_dups,
      unsigned long       n_orig_plants,
      struct PLANT_RECORD *orig_plants_ptr );

   float gauss_dev();
   float uniform_0_1();

   void fill_in_missing_tree_expf(
      unsigned long       *return_code, 
      double              fixed_plot_radius,
      double              min_dbh,
      double              baf,
      unsigned long       n_plants,
      struct PLANT_RECORD *plants_ptr );


/****************************************************************************/
/* functions in species.c                                                   */
/****************************************************************************/

   int compare_species_by_idx(
      const void *ptr1,
      const void *ptr2 );

   int compare_species_by_sp_code(
      const void *ptr1,
      const void *ptr2 );

   struct SPECIES_RECORD   *get_species_entry_from_code(
      unsigned long           n_records,
      struct SPECIES_RECORD   *species_ptr,
      const char              *sp_code );

   struct SPECIES_RECORD *read_species_file(
      unsigned long   *return_code,
      const char      *filename, 
      unsigned long   *n_species );

   void write_species_file(    
      unsigned long           *return_code,
      const char              *filename, 
      unsigned long           n_records,
      struct SPECIES_RECORD   *species_ptr );



/****************************************************************************/
/* functions in stats.c                                                     */
/****************************************************************************/
   void fill_in_whc_and_precip( 
      unsigned long           *return_code,
      double					dfsi,
      double					ppsi,
      double					*whc,
      double					*precip);

   void calc_plot_stats_2( 
      unsigned long           *return_code,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_points,
      struct PLOT_RECORD      *plots_ptr );

   void calc_values_in_taller_2( 
      unsigned long                       *return_code,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long                       n_coeffs,   
      struct COEFFS_RECORD       *coeffs_ptr,
      unsigned long                       n_plants,
      struct PLANT_RECORD         *plants_ptr,
      unsigned long                       n_points,
      struct PLOT_RECORD         *plots_ptr );


   struct SUMMARY_RECORD *build_species_summaries( 
      unsigned long           *return_code,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           *n_sp_in_sample );

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
      struct SUMMARY_RECORD   *summaries_ptr );

/* this function will return a pointer to a structure that has          */
/* the code, all the entries in teh summaries_ptr should have a unique  */
/* code. It could be species codes or functional species code or        */
/* something else                                                       */
   struct SUMMARY_RECORD *get_summary_from_code(
      unsigned  long          n_summaries,
      struct SUMMARY_RECORD   *summaries_ptr,
      unsigned long           code );

   void update_total_summaries( 
      unsigned long           *return_code,
      unsigned long           n_points,
      unsigned long           n_plants,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr,
      struct PLANT_RECORD     *plants_ptr,
      struct SUMMARY_RECORD   *sum_ptr );


   struct SUMMARY_RECORD *build_fsp_summaries( 
      unsigned long           *return_code,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           *n_fsp_in_sample );

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
      struct SUMMARY_RECORD   *summaries_ptr );

   void __stdcall calc_max_sdi(
      unsigned long   *return_code,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      unsigned long           n_points,
      double                  *sdimx  );

   void calc_sites_from_awi( 
      double  awi, 
      double  *df_site, 
      double  *pp_site );

   void calc_sites_from_whc(
      unsigned long   *return_code,
      double          whc,
      double          precip,
      double          *df_site,
      double          *pp_site );

   void get_type_count( 
      unsigned long           *return_code,
      unsigned long           type,
      unsigned long           *type_count,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr );

void set_in_taller_attribs(
    struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *c_ptr,
    struct PLOT_RECORD      *plot_ptr );

void get_in_taller_attribs(
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    double                  *bait,
    double                  *cait );

void set_in_larger_attribs(
   struct PLANT_RECORD     *plant_ptr,
    struct COEFFS_RECORD    *c_ptr,
    struct PLOT_RECORD      *plot_ptr );

void get_in_larger_attribs(
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    double                  *bal);


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
	  double                  baf);

/* this function is deprecated, and calls the previous function */
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
	  double                  baf);





/****************************************************************************/
/* functions in thin.c                                                      */
/****************************************************************************/
   void __stdcall thin_plot( 
      unsigned long           *return_code,
      unsigned long           n_plants,
      struct PLANT_RECORD     *plants_ptr,
      struct PLOT_RECORD      *plot_ptr,
      unsigned long           n_species,
      struct SPECIES_RECORD   *species_ptr,
      unsigned long           n_coeffs,
      struct COEFFS_RECORD    *coeffs_ptr,
      unsigned long           thin_species_idx, 
      int                     thin_guide,
      double                  target, 
      double                  *plants_removed, 
      double                  *ba_removed );


/****************************************************************************/
/* functions in qaqc_mgt.c                                                  */
/****************************************************************************/
   void write_test_file(
        FILE                    *fp,
        unsigned long           sim_variant,
        long                    cycle,

        unsigned long           n_species,
        struct SPECIES_RECORD   *species_ptr,

        unsigned long           n_points,
        struct PLOT_RECORD      *plots_ptr,

        unsigned long           n_plants,
        struct PLANT_RECORD     *plants_ptr,
        
        unsigned long           n_years_after_planting );

   void write_cips_test_file(
    FILE                    *fp,
    long                    cycle,

    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,

    unsigned long           n_points,
    struct PLOT_RECORD      *plots_ptr,

    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    
    unsigned long           n_years_after_planting );

   


#ifdef __cplusplus
}
#endif


#endif /* __CONIFERS_H__ */




