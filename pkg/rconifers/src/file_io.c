/****************************************************************************/
/*                                                                          */
/*  file_io.c                                                               */
/*  functions used to predict values for the CONIFERS growth model          */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                          Revision History                                */
/*                                                                          */
/*  Number  Date        Who     Revision Notes                              */
/****************************************************************************/
/*  MOD000  Sept12,1999 JDH     Created file and added header information   */
/*  MOD001  Sept14,1999 JDH     Removed #include "model.h" and incorped it  */
/*                              into "usfs_r5.h"                            */
/*  MOD002  Sept21,1999 JDH     moved  import_plot_array_from_file() and    */
/*                              write_plot_file() from plot.c to file_io.c  */
/*  MOD003  Oct 18,1999 JDH     added aguments to write_cactos_file()       */
/*                              and write_cactos_ingrowth_file()            */
/*  MOD004  Nov 08,1999 JDH     added plot->lat and plot->long to file i/o  */
/*  MOD005  Nov 13,1999 JDH     added combine_plant_lists()                 */
/*  MOD006  Nov 13,1999 JDH     fixed bug in file i/o functions for invalid */
/*                              problems                                    */
/*  MOD007  Nov 14,1999 JDH     fixed initialization bug in                 */
/*                              read_plant_summary_file()                   */
/*  MOD008  Jan 09,2000 JDH     fixed return_code in                        */
/*                              read_plant_summary_file()                   */
/*  MOD009  Jan 22,2000 JDH     added read_config_record() and              */
/*                              write_config_record()                       */
/*  MOD010  Jan 22,2000 JDH     added copy_sample()                         */
/*  MOD011  Jan 22,2000 JDH     rewrote write_cactos_file()                 */
/*                              and write_fvs_file()                        */
/*  MOD012  Jan 25,2000 MWR     added changes for fvs output                */
/*  MOD013  Jan 25,2000 MWR     added changes for writing a conifers        */
/*                              txt file                                    */
/*  MOD014  Jan 29,2000 mwr     changed write format for conifers file      */
/*                                   added more sig figures for a couple    */
/*  MOD015  Feb 03,2000 MWR     added vars to the conifers file format      */
/*  MOD016  Feb 10,2000 MWR     added read systum archive function          */
/*  MOD017  Feb 11,2000 MWR     added convert systum tree species           */
/*  MOD018  Feb 11,2000 MWR     added convert systum shrub species          */
/*  MOD019  Feb 25,2000 MWR     called the fill_in_whc_precip               */
/*  MOD020  Feb 25,2000 JDH     added write_organon_inp_file()              */
/*  MOD021  Feb 26,2000 JDH     updated write_cactos_ingrowth_file()        */
/*                              to only export tree records                 */
/*  MOD022  Aug 10,2000 JDH     updated read_config_record() and            */
/*                              write_config_record                         */
/*  MOD023  Oct 26,2000 JDH     updated write_organon_inp() to not export   */
/*                              tree records where the expf < 0.0001        */
/*  MOD024  Oct 03,2000 JDH     updated write_summaries_to_file() to        */
/*                              include h/d6 ratios                         */
/*  MOD025  Oct 11,2001 JDH     changed dbh filter for organon output from  */
/*                              DBH <= 0.0 to DBH < 0.1 since organon needs */
/*                              trees > 4.5 and dbh => .1                   */
/*  MOD026  Dec 03,2001 JDH     fixed save as ((double)elevation in sample  */
/*  MOD027  Dec 10,2001 MWR     fixed argument list in function             */
/*  MOD028  Jul 24,2003 MWR    error checking on cr                        */
/****************************************************************************/

/****************************************************************************/
/*  functions list in file_io.c                                             */
/*      1. struct PLANT_RECORD  *read_systum1_file()                        */
/*      2. void                 write_systum1_file( )                       */
/*      3. void                 dump_plots_to_file( )                       */
/*      4. void                 dump_plants_to_file( )                      */
/*      5. void                 write_summaries_to_file( )                  */
/*      6. struct PLOT_RECORD   *read_plots_from_file( )                    */
/*      7. void                 write_plots_to_text_file( )                 */
/*      8. struct SAMPLE_RECORD *read_sample_from_file( )                   */
/*      9. void                 write_sample_to_file( )                     */
/*     10. void                 free_sample_data( )                         */
/*     11. void                 write_organon_file( )                       */ 
/*     12. void                 write_cactos_file(  )                       */
/*     13. void                 write_cactos_ingrowth_file( )               */
/*     14. void                 write_fvs_file( )                           */
/*     15. static int           get_fvs_cr( )                               */
/*     16. struct PLANT_RECORD  *read_usfs_r5_ref_file( )                   */
/*     17. struct PLANT_RECORD  *read_plant_summary_file( )                 */
/*     18. struct PLANT_RECORD  *combine_plant_lists( )                     */
/*     19. void                 read_config_record( )                       */
/*     20. void                 write_config_record( )                      */
/*     21. void                 backup sample( )                            */
/*     22. void                 return_sample( )                            */
/*     23. struct SAMPLE_RECORD read_systum_archive()                       */
/*     24. convert_systum_tree_species                                      */
/*     25. convert_systum_shrub_species                                     */
/*     26. write_organon_inp                                                */
/*                                                                          */
/****************************************************************************/

/* 	$Id: file_io.c 853 2012-01-24 02:05:20Z hamannj $	 */

/* #ifndef lint */
/* static char vcid[] = "$Id: file_io.c 853 2012-01-24 02:05:20Z hamannj $"; */
/* #endif /\* lint *\/ */


//#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"

/* local functions */
/* static int compare_plants_by_plot_and_sp(  */
/*     const void *ptr1,  */
/*     const void *ptr2 ); */

static int get_fvs_cr( double cr );

static unsigned long get_org_sp_group( 
    unsigned long           *return_code,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    struct PLANT_RECORD     *plant_ptr,
    unsigned long           *is_big6 );

/* MOD052 */
static void convert_systum_tree_species( 
      unsigned long           *return_code,
      double                  d_code,
      char                    *spcode);

/* MOD052 */ 
static void convert_systum_shrub_species( 
      unsigned long           *return_code,
      double                  d_code,
      char                    *spcode);

/* moved from species.c */

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


/****************************************************************************/
/* sorting functions for species records                                    */
/****************************************************************************/

int compare_species_by_idx(
    const void *ptr1,
    const void *ptr2 )
{
    struct SPECIES_RECORD   *sp1_ptr;
    struct SPECIES_RECORD   *sp2_ptr;

    sp1_ptr = (struct SPECIES_RECORD*)ptr1;
    sp2_ptr = (struct SPECIES_RECORD*)ptr2;

    if( sp1_ptr->idx > sp2_ptr->idx )
    {
        return 1;
    }
    if( sp1_ptr->idx < sp2_ptr->idx )
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

int compare_species_by_sp_code(
    const void *ptr1,
    const void *ptr2 )
{
    struct SPECIES_RECORD   *sp1_ptr;
    struct SPECIES_RECORD   *sp2_ptr;

    sp1_ptr = (struct SPECIES_RECORD*)ptr1;
    sp2_ptr = (struct SPECIES_RECORD*)ptr2;

    /* stricmp isn't ansi!                                      */
    /* return stricmp( sp1_ptr->sp_code, sp2_ptr->sp_code );    */
    return strcmp( sp1_ptr->sp_code, sp2_ptr->sp_code );

}



struct SPECIES_RECORD   *get_species_entry_from_code(
    unsigned long           n_records,
    struct SPECIES_RECORD   *species_ptr,
    const char              *sp_code )
{

   struct SPECIES_RECORD   *entry;
   struct SPECIES_RECORD   key;
   
   /* this is a runtime check */
   if( sp_code == NULL )
   {
      return NULL;
   }
   
   memset( &key, 0, sizeof( struct SPECIES_RECORD ) );
   strcpy( key.sp_code, sp_code );
   
   /* find the matching species record  */
   entry = (struct SPECIES_RECORD*)bsearch(
      &key,
      species_ptr,
      n_records,
      sizeof( struct SPECIES_RECORD ),
      compare_species_by_sp_code );
   
    return entry;
    
}






/****************************************************************************/
/* 1. read_systum1_file                                                     */
/****************************************************************************/
void read_systum1_file( 
    unsigned long           *return_code,
    const char              *filename, 
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           *n_points,
    struct PLOT_RECORD      **plots_ptr,
    unsigned long           *n_plants,
    struct PLANT_RECORD     **plants_ptr,
	unsigned long			*age )	
{

    char                line_buffer[256];
    unsigned long       i;
	long                temp_plot;
	char                temp_species[SP_LENGTH];
	double	            temp_dbh;
	double	            temp_tht;
	double	            temp_cr;
	double	            temp_expf;
    unsigned long       temp_uc;
    //struct PLANT_RECORD *tree_ptr;
	FILE	            *fp;

    char                    temp_sp_code[16];
    struct SPECIES_RECORD   *s_ptr;

    /* open the file */
    if( ( fp = fopen( filename, "rt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }


    /* sort the species codes based on sp_code */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_sp_code );


	if( *plots_ptr )
	{
		free( *plots_ptr );
		*plots_ptr = NULL;
	}

	if( *plants_ptr )
	{
		free( *plants_ptr );
		*plants_ptr = NULL;
	}


    /*  count the number of lines in the file   */
    /*  1 PP -1 2.2 -1 100                      */    
	i = 0;
	while( fgets( line_buffer, sizeof( line_buffer ), fp ) != NULL )
	{
		if( sscanf( line_buffer,
		        "%ld %s %lf %lf %lf %lf %ld", 
				&temp_plot,
				temp_species,
				&temp_dbh,
				&temp_tht,
				&temp_cr,
				&temp_expf,
				&temp_uc ) < 6 )
		{
			fclose( fp );
            *return_code = INVALID_OPTION;
            (*n_plants) = 0;
            //return NULL;
            return;
			//continue;
		}

		i++;
	}

    (*n_plants) = i;
    *plants_ptr = (struct PLANT_RECORD *)calloc( 
        (*n_plants), sizeof( struct PLANT_RECORD ) );

    rewind( fp );
    
    i = 0;
	/* read each product for the rest of the file   */
	while( fgets( line_buffer, sizeof( line_buffer ), fp ) != NULL )
	{
        memset( &(*plants_ptr)[i], 0, sizeof( struct PLANT_RECORD ) );

        if( sscanf( line_buffer, 
		        "%ld %s %lf %lf %lf %lf %ld", 
				&(*plants_ptr)[i].plot,
				//tree_ptr[i].sp_code,
                temp_sp_code,
				&(*plants_ptr)[i].dbh,
				&(*plants_ptr)[i].tht,
				&(*plants_ptr)[i].cr,
				&(*plants_ptr)[i].expf,
				&(*plants_ptr)[i].user_code ) < 6 )
        {
            continue;/* should never get here because of trap on first read above */
        }


        /* convert the temporary species code to the internal functional species code */
        //tree_ptr[i].sp_idx = get_sp_idx( temp_sp_code );
        s_ptr = get_species_entry_from_code(    n_species,
                                                species_ptr, 
                                                temp_sp_code );

        if( !s_ptr )
        {
            *return_code = INVALID_OPTION;
            return;
        }

        (*plants_ptr)[i].sp_idx = s_ptr->idx;

		/* do some crude fixing here    */
		if( (*plants_ptr)[i].dbh < 0.0 )
        {
            (*plants_ptr)[i].dbh = 0.0;
        }
		if( (*plants_ptr)[i].tht < 0.0 )
        {
            (*plants_ptr)[i].tht = 0.0;
        }
		if( (*plants_ptr)[i].cr < 0.0 )
        {
            (*plants_ptr)[i].cr = 0.0;
        }
		/*  MOD028 */
		if( (*plants_ptr)[i].cr > 1.0)
		{
			(*plants_ptr)[i].cr =1.0;
		}

		if( (*plants_ptr)[i].expf < 0.0 )
        {
            (*plants_ptr)[i].expf = 0.0;
        }
        
        (*plants_ptr)[i].n_stems  = 1;
        
        //tree_ptr[i]->basal_area = tree_ptr[i]->dbh * tree_ptr[i]->dbh * FC_I;
        (*plants_ptr)[i].basal_area = (*plants_ptr)[i].dbh * (*plants_ptr)[i].dbh * FC_I;

        i++;    
    }

    fclose( fp );
    
    *return_code = CONIFERS_SUCCESS;

    /* then resort the species codes based on idx */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_idx );

	//return tree_ptr;

}




/****************************************************************************/
/* 2. write_systum1_file                                                    */
/****************************************************************************/
void __stdcall write_systum1_file( 
    unsigned long       *return_code,
    const char          *filename, 
    unsigned long       n_trees,
    struct PLANT_RECORD *tree_ptr,
    unsigned long       n_species,
    struct SPECIES_RECORD *species_ptr )
{

	FILE	            *fp;
    struct PLANT_RECORD *T;
    unsigned long       i;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    T = &tree_ptr[0];
    for( i = 0; i < n_trees; i++, T++ )
    {
        fprintf( fp, 
            "%ld %s %lf %lf %lf %lf %ld\n", 
            T->plot,
            //T->sp_code,
            species_ptr[T->sp_idx].sp_code,
            T->dbh,
            T->tht,
            T->cr,
            T->expf,
            T->user_code );
    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;
}


/****************************************************************************/
/* 3. dump_plots_to_file                                                    */
/****************************************************************************/
void __stdcall dump_plots_to_file( 
    unsigned long       *return_code,
    const char          *filename,
    unsigned long       n_points,
    struct PLOT_RECORD  *plots_ptr )
{

	FILE	            *fp;
    struct PLOT_RECORD  *plot_ptr;
    unsigned long       i;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

        fprintf( fp, 
//            "%4ld %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf\n", 
            "plot, water_capacity, site_30, shrub_pct_cover, shrub_mean_height, basal_area, bh_expf, qmd, sdi, ba_c, ba_h, expf\n" );

    plot_ptr = &plots_ptr[0];
    for( i = 0; i < n_points; i++, plot_ptr++ )
    {
        /* MOD004 */
        fprintf( fp, 
//            "%4ld %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf\n", 
            "%4ld, %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf\n", 
            plot_ptr->plot,          
            plot_ptr->water_capacity, 
            plot_ptr->site_30,
            plot_ptr->shrub_pct_cover,  
            plot_ptr->shrub_mean_height,
            plot_ptr->basal_area,
            plot_ptr->bh_expf,               
            plot_ptr->qmd,              
            plot_ptr->sdi,
            plot_ptr->ba_c,
            plot_ptr->ba_h,
            plot_ptr->expf

			);

    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

}



/****************************************************************************/
/*  4. dump_plants_to_file                                                  */
/****************************************************************************/
void __stdcall dump_plants_to_file( 
    unsigned long       *return_code,
    const char          *filename,
//    FILE                *fp,
    unsigned long       n_plants,
    struct PLANT_RECORD *plants_ptr,
    unsigned long       n_species,
    struct SPECIES_RECORD *species_ptr )
{

	FILE	            *fp;
    struct PLANT_RECORD *plant_ptr;
    unsigned long       i;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
       *return_code = INVALID_OPTION;
        return;
    }

        fprintf( fp, 
            //"%ld %ld %s %lf %lf %lf %lf %lf %lf %ld %lf %lf %lf %ld %lf %lf %lf %lf %lf %lf %ld\n", 
            "plot, plant, sp_code, d6, d6_area, dbh, basal_area, tht, cr, n_stems, expf, crown_width, crown_area, user_code, d6_growth, dbh_growth, tht_growth, cr_growth, cw_growth, expf_change, errors\n" );


    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        fprintf( fp, 
            //"%ld %ld %s %lf %lf %lf %lf %lf %lf %ld %lf %lf %lf %ld %lf %lf %lf %lf %lf %lf %ld\n", 
            "%ld, %ld, %s, %lf, %lf, %lf, %lf, %lf, %lf, %ld, %lf, %lf, %lf, %ld, %lf, %lf, %lf, %lf, %lf, %lf, %ld\n", 
            plant_ptr->plot,            
            plant_ptr->plant,           

            //plant_ptr->sp_code,         
            species_ptr[plant_ptr->sp_idx].sp_code,

            plant_ptr->d6,              
            plant_ptr->d6_area,         
            plant_ptr->dbh,             
            plant_ptr->basal_area,      
            plant_ptr->tht,             
            plant_ptr->cr,			    
            plant_ptr->n_stems,         
            plant_ptr->expf,            
            plant_ptr->crown_width,     
            plant_ptr->crown_area,      
            plant_ptr->user_code,       
            plant_ptr->d6_growth,       
            plant_ptr->dbh_growth,      
            plant_ptr->tht_growth,      
            plant_ptr->cr_growth,       
            plant_ptr->cw_growth,       
            plant_ptr->expf_change,     
            plant_ptr->errors );
            //plant_ptr->bait_c,          
            //plant_ptr->bait_h,          
            //plant_ptr->bait_s );        

    }

//    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

}


/****************************************************************************/
/*  5. write_summaries_to_file                                              */
/****************************************************************************/
void __stdcall write_summaries_to_file( 
    unsigned long           *return_code,
    const char              *filename,
    unsigned long           n_records,
    struct SUMMARY_RECORD   *summaries_ptr,
    unsigned long			n_species,
    struct SPECIES_RECORD	*species_ptr,
    unsigned long           age)
{

	FILE	                *fp;
    struct SUMMARY_RECORD   *sum_ptr;
    unsigned long           i;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    fprintf(    fp, 
                "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", 
                " age",
                "sp_code",
                "min_dbh",               
                "qmd",                   
                "max_dbh",               
                "basal_area",            
                "min_height",            
                "mean_height",           
                "max_height",            

                "min_hd6_ratio",
                "mean_hd6_ratio",
                "max_hd6_ratio",

                "expf",                  
                "cr",
                "crown_area",            
                "pct_cover",             
                "sdi" );

    sum_ptr = &summaries_ptr[0];
    for( i = 0; i < n_records; i++, sum_ptr++ )
    {
        fprintf( fp, 
            "%ld,%s,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
            age,                    /*  time or age for the summaries   */
            species_ptr[sum_ptr->code].sp_code,
            sum_ptr->min_dbh,       /*  minumim dbh                     */
            sum_ptr->qmd,           /*  quad mean diameter              */
            sum_ptr->max_dbh,       /*  max dbh in the list             */
            sum_ptr->basal_area,    /*  basal area for the code         */
            sum_ptr->min_height,    /*  min total height                */
            sum_ptr->mean_height,   /*  average total height            */
            sum_ptr->max_height,    /*  max total height                */
            sum_ptr->min_hd6_ratio, /*  min total height                */
            sum_ptr->mean_hd6_ratio,/*  average total height            */
            sum_ptr->max_hd6_ratio, /*  max total height                */
            sum_ptr->expf,          /*                                  */
            sum_ptr->cr,            /*                                  */
            sum_ptr->crown_area,    /*                                  */
            sum_ptr->pct_cover,     /*                                  */
            sum_ptr->sdi );			/*                                  */
    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

}



/****************************************************************************/
/* read_plots_from_file												*/
/* file mimics the read_plots_from_file function							*/
/****************************************************************************/
void read_plots_from_file( 
	unsigned long			*return_code,
    const char				*filename, 
    unsigned long			*n_plots,
	struct PLOT_RECORD		**plots_ptr )
{

    char    line_buffer[256];
	long    i;
	FILE	*fp;
    //struct PLOT_RECORD *plot_ptr;

    if( ( fp = fopen( filename, "rt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

	/* eat the header line */
	fgets( line_buffer, sizeof( line_buffer ), fp );

	i = 0;
	while( fgets( line_buffer, sizeof( line_buffer ), fp ) != NULL )
	{
		i++;
	}

    (*n_plots) = i;
    *plots_ptr = (struct PLOT_RECORD *)calloc( 
        (*n_plots), sizeof( struct PLOT_RECORD ) );

    rewind( fp );

	/* eat the header line */
	fgets( line_buffer, sizeof( line_buffer ), fp );

    i = 0;
	/*  read each product for the rest of the file */
	while( fgets( line_buffer, sizeof( line_buffer ), fp ) != NULL )
	{
        memset( &(*plots_ptr)[i], 0, sizeof( struct PLOT_RECORD ) );

        /* MOD004 */
        if( sscanf( line_buffer, 
		        "%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf"
				"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
				"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",

				&(*plots_ptr)[i].plot,
				&(*plots_ptr)[i].latitude,
				&(*plots_ptr)[i].longitude,
				&(*plots_ptr)[i].elevation,
				&(*plots_ptr)[i].slope,
				&(*plots_ptr)[i].aspect,
				&(*plots_ptr)[i].water_capacity,
				&(*plots_ptr)[i].mean_annual_precip,
				&(*plots_ptr)[i].site_30,
				
				/* added for swohybrid */
				&(*plots_ptr)[i].growing_season_precip, 

				/* monthly average temperatures */
				&(*plots_ptr)[i].mean_monthly_temp[0],
				&(*plots_ptr)[i].mean_monthly_temp[1],
				&(*plots_ptr)[i].mean_monthly_temp[2],
				&(*plots_ptr)[i].mean_monthly_temp[3],
				&(*plots_ptr)[i].mean_monthly_temp[4],
				&(*plots_ptr)[i].mean_monthly_temp[5],
				&(*plots_ptr)[i].mean_monthly_temp[6],
				&(*plots_ptr)[i].mean_monthly_temp[7],
				&(*plots_ptr)[i].mean_monthly_temp[8],
				&(*plots_ptr)[i].mean_monthly_temp[9],
				&(*plots_ptr)[i].mean_monthly_temp[10],
				&(*plots_ptr)[i].mean_monthly_temp[11],

				/* solar radiation */
				&(*plots_ptr)[i].solar_radiation[0],
				&(*plots_ptr)[i].solar_radiation[1],
				&(*plots_ptr)[i].solar_radiation[2],
				&(*plots_ptr)[i].solar_radiation[3],
				&(*plots_ptr)[i].solar_radiation[4],
				&(*plots_ptr)[i].solar_radiation[5],
				&(*plots_ptr)[i].solar_radiation[6],
				&(*plots_ptr)[i].solar_radiation[7],
				&(*plots_ptr)[i].solar_radiation[8],
				&(*plots_ptr)[i].solar_radiation[9],
				&(*plots_ptr)[i].solar_radiation[10],
				&(*plots_ptr)[i].solar_radiation[11]
				
				) < 9 )
				//&(*plots_ptr)[i].site_30 ) < 7 )
        {
            continue;/* error trap here */
        }

        i++;    
    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

}


/****************************************************************************/
/*  7. write_plot_file                                                      */
/****************************************************************************/
/* MODXXX */
/* this function is deprecated since it's assumed the user will be using	*/
/* the rconifers library from R												*/
void __stdcall write_plots_to_text_file( 
    unsigned long       *return_code,
    const char          *filename, 
    unsigned long       n_records,
    struct PLOT_RECORD  *plots_ptr )
{

	FILE	            *fp;
    struct PLOT_RECORD  *plot_ptr;
    unsigned long       i;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    plot_ptr = &plots_ptr[0];
    for( i = 0; i < n_records; i++, plot_ptr++ )
    {
        /* MOD004 */
        fprintf( fp, 
            "%ld %lf %lf %lf %lf %lf %lf %lf %lf\n", 
            plot_ptr->plot,
            plot_ptr->latitude,
            plot_ptr->longitude,
            plot_ptr->elevation,
            plot_ptr->slope,
            plot_ptr->aspect,
            plot_ptr->water_capacity,
            plot_ptr->mean_annual_precip,
            plot_ptr->site_30 );
    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;
}


/****************************************************************************/
/* this function reads a cleaned up sample file                             */
/* from the disk, and allocates the SAMPLE_RECORD                           */
/* arrays and fills them in                                                 */
/*   this is essentially a conifers file                                    */
/****************************************************************************/
/*   8. read_sample_from_file                                               */
/****************************************************************************/
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
    unsigned long           *variant)

{
	FILE	                *fp;
    struct PLOT_RECORD      *plot_ptr = NULL;
    struct PLANT_RECORD     *plant_ptr = NULL;
    //struct SAMPLE_RECORD    *sample_ptr;
    unsigned long           i;
    int                     n_args;

    char    line_buffer[256];

    char                    temp_sp_code[16];
    struct SPECIES_RECORD   *s_ptr;

    if( ( fp = fopen( filename, "rt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    //sample_ptr = (struct SAMPLE_RECORD*)calloc( 1, sizeof( struct SAMPLE_RECORD ) );

    //if( sample_ptr == NULL )
    //{
    //    *return_code = CONIFERS_ERROR;
    //    return NULL;
   // }

    /* clean up any previous pointers */
    /* and data */
    if( *plots_ptr )
    {
        free( *plots_ptr );
        *plots_ptr = NULL;
    }

    if( *plants_ptr )
    {
        free( *plants_ptr );
        *plants_ptr = NULL;
    }

    
    /*  MOD012   */
    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, "%60s", &sample_ptr->forest );
    //strncpy( sample_ptr->forest, line_buffer, strlen( line_buffer ) - 1 );

    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, "%60s", &sample_ptr->subunit);
    //strncpy( sample_ptr->subunit, line_buffer, strlen( line_buffer ) - 1 );

    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, "%60s", &sample_ptr->stand_name);
    //strncpy( sample_ptr->stand_name, line_buffer, strlen( line_buffer ) - 1 );

    /* this is lame */
    /* you have to take all the cr/lf's out */
    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, "%60s", &sample_ptr->legal);
    //strncpy( sample_ptr->legal, line_buffer, strlen( line_buffer ) - 1 );


    /*  MOD015 */
    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, "%ld", &sample_ptr->elevation);

    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, "%lf",  &sample_ptr->acreage);

    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, "%ld",  &sample_ptr->age);
    n_args = sscanf( line_buffer, "%ld %ld",  age, variant);
    if(n_args != 2)
    {
        *variant=CONIFERS_SWO;
    }
    if(*variant < CONIFERS_SWO || *variant > CONIFERS_SMC)
    {
        *return_code = INVALID_VARIANT;
        return;   
    }

    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, "%ld", &sample_ptr->current_year);

/*
    fgets( line_buffer, sizeof( line_buffer ), fp );    
    sscanf( line_buffer, "%12s%ld%ld", &sample_ptr->month,
                            &sample_ptr->day,
                            &sample_ptr->year);
*/

    fgets( line_buffer, sizeof( line_buffer ), fp );    
    //sscanf( line_buffer, 
    //        "%ld%ld%ld", 
    //        &sample_ptr->sampled_month,
    //        &sample_ptr->sampled_day,
    //        &sample_ptr->sampled_year);

    /* read the first line of the file to get the number of plots */    
	fgets( line_buffer, sizeof( line_buffer ), fp );
    //sscanf( line_buffer, "%ld", &sample_ptr->n_points );
    sscanf( line_buffer, "%ld", n_points );

/*
    sample_ptr->plots_ptr = (struct PLOT_RECORD*)calloc( 
        sample_ptr->n_points, sizeof( struct PLOT_RECORD ) );
*/



    *plots_ptr = (struct PLOT_RECORD*)calloc( 
		*n_points, sizeof( struct PLOT_RECORD ) );

    //plot_ptr = &sample_ptr->plots_ptr[0];
    plot_ptr = &(*plots_ptr)[0];
    //for( i = 0; i < sample_ptr->n_points; i++, plot_ptr++ )
    for( i = 0; i < (*n_points); i++, plot_ptr++ )
    {
	    fgets( line_buffer, sizeof( line_buffer ), fp );
        
        if( sscanf( line_buffer, 
		        "%ld %lf %lf %lf %lf %lf %lf %lf %lf", 
				&plot_ptr->plot,
                &plot_ptr->latitude,
                &plot_ptr->longitude,
				&plot_ptr->elevation,
				&plot_ptr->slope,
				&plot_ptr->aspect,
				&plot_ptr->water_capacity,
				&plot_ptr->mean_annual_precip,
				&plot_ptr->site_30 ) < 6 )
        {
            *return_code = CONIFERS_ERROR;


			if( *plots_ptr )
			{
				free( *plots_ptr );
                *plots_ptr = NULL;
			}

            //if( sample_ptr->plots_ptr )
            //{
            //    free( sample_ptr->plots_ptr  );
            //}
            //if( sample_ptr )
            //{
            //    free( sample_ptr );
            //}

            fclose( fp );

            //return NULL;
            return;
        }

    }        


    /* sort the species codes based on sp_code */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_sp_code );


    /* now import the plants in the file */    
	fgets( line_buffer, sizeof( line_buffer ), fp );
    //sscanf( line_buffer, "%ld", &sample_ptr->n_plants );
    sscanf( line_buffer, "%ld", n_plants );

    //sample_ptr->plants_ptr = (struct PLANT_RECORD*)calloc( 
    //    sample_ptr->n_plants, sizeof( struct PLANT_RECORD ) );

    *plants_ptr = (struct PLANT_RECORD*)calloc( 
        *n_plants, sizeof( struct PLANT_RECORD ) );


    //plant_ptr = &sample_ptr->plants_ptr[0];
    //for( i = 0; i < sample_ptr->n_plants; i++, plant_ptr++ )


    plant_ptr = &(*plants_ptr)[0];
    for( i = 0; i < (*n_plants); i++, plant_ptr++ )
    {

      n_args = 0;
      memset( line_buffer, 0, sizeof( line_buffer ) );
      memset( temp_sp_code, 0, sizeof( temp_sp_code ) );

	    fgets( line_buffer, sizeof( line_buffer ), fp );
        
        n_args = sscanf( line_buffer, 
            "%ld %ld %s %lf %lf %lf %lf %lf %lf %ld %lf %lf %lf %ld", 
            &plant_ptr->plot,           
            &plant_ptr->plant,          
            temp_sp_code,
            &plant_ptr->d6,             
            &plant_ptr->d6_area,        
            &plant_ptr->dbh,            
            &plant_ptr->basal_area,     
            &plant_ptr->tht,            
            &plant_ptr->cr,			    
            &plant_ptr->n_stems,        
            &plant_ptr->expf,           
            &plant_ptr->crown_width,    
            &plant_ptr->crown_area,     
            &plant_ptr->user_code  );

        /* convert the temporary species code to the internal functional species code */
        //tree_ptr[i].sp_idx = get_sp_idx( temp_sp_code );
        s_ptr = NULL;
        s_ptr = get_species_entry_from_code(    n_species,
                                                species_ptr, 
                                                temp_sp_code );


        /* if the species is invalid for the current species list */
        /* clean up the sample and return null, setting the       */
        /* rc to invalid file                                     */
        if( !s_ptr || n_args < 14 )
        {

            *return_code = CONIFERS_ERROR;

			if( *plots_ptr )
			{
				free( *plots_ptr );
                *plots_ptr = NULL;

			}

			if( *plants_ptr )
			{
				free( *plants_ptr );
                *plants_ptr = NULL;
			}

            //if( sample_ptr->plots_ptr )
            //{
            //    free( sample_ptr->plots_ptr  );
            //}

            //if( sample_ptr->plants_ptr )
            //{
            //    free( sample_ptr->plants_ptr  );
            //}

            //if( sample_ptr )
            //{
            //    free( sample_ptr );
            //}


            fclose( fp );

            //return NULL;
			return;

        }

        plant_ptr->sp_idx = s_ptr->idx;

    }

    fclose( fp );

    *return_code = CONIFERS_SUCCESS;

    /* sort the species codes based on idx */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_idx );

    //return sample_ptr;
	return;

}


/****************************************************************************/
/*   9.  write_sample_to_file                                               */
/****************************************************************************/
void __stdcall write_sample_to_file( 
    unsigned long           *return_code,
    const char              *filename, 
    //struct SAMPLE_RECORD    *sample_ptr,
    unsigned long			age,
    unsigned long           variant,
    unsigned long           n_points,
    struct PLOT_RECORD      *plots_ptr,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr )
{

	FILE	            *fp;
    struct PLOT_RECORD  *plot_ptr;
    struct PLANT_RECORD *plant_ptr;
    unsigned long       i;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }
    
    /* MOD012  */
    /* MOD015  */
    /*  print out the header stuff  */
    /* MOD026 cast this to a double */
    //fprintf( fp, "%s\n",    sample_ptr->forest);
    //fprintf( fp, "%s\n",    sample_ptr->subunit);
    //fprintf( fp, "%s\n",    sample_ptr->stand_name);
    //fprintf( fp, "%s\n",    sample_ptr->legal);
	//fprintf( fp, "%lf\n",   (double)sample_ptr->elevation);     
    //fprintf( fp, "%lf\n",    sample_ptr->acreage);
    //fprintf( fp, "%ld\n",           sample_ptr->current_year);
    
    fprintf( fp, "%s\n",    "## (forest) deprecated." );
    fprintf( fp, "%s\n",    "## (subunit) deprecated." );
    fprintf( fp, "%s\n",    "## (stand name) deprecated." );
    fprintf( fp, "%s\n",    "## (legal desription) deprecated." );
    fprintf( fp, "%s\n",    "## (mean stand elevation) deprecated." );
    fprintf( fp, "%s\n",    "## (acres) deprecated." );
    fprintf( fp, "%ld %ld\n",       age, variant);
    fprintf( fp, "%s\n",    "## (current year) deprecated." );
    fprintf( fp, "%s\n",    "## (sampled mm/dd/yyyy) deprecated." );


    //fprintf( fp, "%ld\n",           sample_ptr->current_year);
    //fprintf( fp, "%s%4ld%6ld\n",    sample_ptr->month,
    //                                sample_ptr->day,
    //                                sample_ptr->year);
    //fprintf( fp, "%ld%4ld%6ld\n",    sample_ptr->sampled_month,
    //                                sample_ptr->sampled_day,
    //                                sample_ptr->sampled_year);

    //fprintf( fp, "%ld\n", sample_ptr->n_points );
    //plot_ptr = &sample_ptr->plots_ptr[0];
    fprintf( fp, "%ld\n", n_points );
    plot_ptr = &plots_ptr[0];

    //for( i = 0; i < sample_ptr->n_points; i++, plot_ptr++ )
    for( i = 0; i < n_points; i++, plot_ptr++ )
    {
        /* MOD004 */
        fprintf( fp, 
            "%5ld %lf %lf %7.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf\n", 
            plot_ptr->plot,
            plot_ptr->latitude,
            plot_ptr->longitude,
            plot_ptr->elevation,
            plot_ptr->slope,
            plot_ptr->aspect,
            plot_ptr->water_capacity,
            plot_ptr->mean_annual_precip,
            plot_ptr->site_30 );
    }


/*   MOD014  changed format statement in fprintf below */
//    fprintf( fp, "%ld\n", sample_ptr->n_plants );
//    plant_ptr = &sample_ptr->plants_ptr[0];
//    for( i = 0; i < sample_ptr->n_plants; i++, plant_ptr++ )

    fprintf( fp, "%ld\n", n_plants );
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        fprintf( fp, 
            "%5ld %5ld %6s %6.2lf %8.4lf %6.2lf %8.4lf %6.2lf %6.3lf %5ld %8.2lf %6.2lf %6.2lf %5ld\n", 
            plant_ptr->plot,            
            plant_ptr->plant,           

            //plant_ptr->sp_code,         
            species_ptr[plant_ptr->sp_idx].sp_code,
            
            
            plant_ptr->d6,              
            plant_ptr->d6_area,         
            plant_ptr->dbh,             
            plant_ptr->basal_area,      
            plant_ptr->tht,             
            plant_ptr->cr,			    
            plant_ptr->n_stems,         
            plant_ptr->expf,            
            plant_ptr->crown_width,     
            plant_ptr->crown_area,      
            plant_ptr->user_code );     
    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

}



/****************************************************************************/
/*  10. free_sample                                                         */
/****************************************************************************/
void free_sample_data( 
    unsigned long           *n_points,
    struct PLOT_RECORD      **plots_ptr,
    unsigned long           *n_plants,
    struct PLANT_RECORD     **plants_ptr )
{

    if( *plots_ptr )
    {
        free( *plots_ptr );
        *n_points     = 0;
        *plots_ptr   = NULL;
    }

    if( *plants_ptr )
    {
        free( *plants_ptr );
        *n_plants    = 0;
        *plants_ptr  = NULL;
    }

}



/****************************************************************************/
/*  11. write_organon_file                                                  */
/****************************************************************************/
void __stdcall write_organon_file( 
    unsigned long           *return_code,
    const char              *filename, 
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    struct COEFFS_RECORD    *coeffs_ptr)
{

	FILE	                *fp;
    struct PLANT_RECORD     *plant_ptr;
    unsigned long           i;
    unsigned long           last_blank;
    unsigned long           last_plot;
    unsigned long           empty_plot;
    unsigned long           empty_plot_num;

    last_blank=0;
    last_plot=0;
    empty_plot=0;
    empty_plot_num=0;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        if(  is_non_stocked(&coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx] ) /* non stocked       */
          || !species_ptr[plant_ptr->sp_idx].organon_sp_code                      /* not organon plant */      
          || !is_tree( &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx])       /* not a tree        */ 
          || plant_ptr->dbh  < 0.1                                                /* too little        */
          || plant_ptr->tht  <= 4.5                                                /* too short         */ 
          || plant_ptr->expf < 0.00001 )                                          /* expf is ~zero     */ 
        {    /* then it should be skipped unless we need to write a placeholder for empty plot         */
             if(last_plot == plant_ptr->plot)
             {
                 continue;
             }
             /* then set a flag                     */
             /* check to see if last was also bogus */
             if(empty_plot && empty_plot_num != plant_ptr->plot)
             {
                    fprintf(fp, "%3ld\n", empty_plot_num);
             }
             empty_plot     = 1;
             empty_plot_num = plant_ptr->plot;
             last_blank     = 1;
             //last_plot  = plant_ptr->plot;
        }
        else 
        {
             /*question: do we need to write out a empty plot in buffer? */
             if(empty_plot && empty_plot_num != plant_ptr->plot)
             {
                    fprintf(fp, "%3ld\n", empty_plot_num);
             }
             fprintf(fp, 
                    "%3ld %3ld %5.1lf %5.1lf %4.2lf %6.2lf %2ld\n",
                    plant_ptr->plot,
                    species_ptr[plant_ptr->sp_idx].organon_sp_code,
			        plant_ptr->dbh,
			        plant_ptr->tht,
			        plant_ptr->cr,
			        plant_ptr->expf,
                    plant_ptr->user_code );

             last_blank = 0;
             last_plot  = plant_ptr->plot;
             empty_plot = 0;
        }
    }

    if(empty_plot)
    {
        fprintf(fp, "%3ld\n", empty_plot_num);
    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

}





/* MOD003 */
/* MOD011 */
/****************************************************************************/
/*  12. write_cactos_file                                                   */
/****************************************************************************/
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
    double                  *site_indecies_ptr,     /* MOD011   */
    unsigned long           n_ages,
    double                  *ages_ptr,              /* MOD011   */
    unsigned long           n_coeffs,               /* MOD011   */
    struct COEFFS_RECORD    *coeffs_ptr )           /* MOD011   */
{

	FILE	                *fp;
    struct PLANT_RECORD     *plant_ptr;
    unsigned long           i;
    unsigned long           total_tree_records;


    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    /* get the count of the tree records in tha plant list */
    total_tree_records = 0;
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        
        if( is_tree( &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx] ) )
        {
            total_tree_records++;
        }

    }
    

    /* export the first line of the file */
	fprintf(    fp, 
/*                 "%-20s%10i\n",  */
                "%-20s%10ld\n",
                sample_id, 
                total_tree_records );

    /* MOD003   */
    /* write the 50 year site index values and ages for the     */
    /* 13 functional species CACTOS can handle                  */
    for( i = 0; i < n_site_indecies; i++ )
    {
	    fprintf(    fp, "%5.0lf.", site_indecies_ptr[i] ); 
    }
    fprintf(    fp, "\n" ); 

    for( i = 0; i < n_ages; i++ )
    {
	    fprintf(    fp, "%5.0lf.", ages_ptr[i] ); 
    }
    fprintf(    fp, "\n" ); 

    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        /* if the plant is not a tree, then don't export it */
        if( !is_tree( &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx] ) )
        {
            continue;
        }

        fprintf(    fp, 
                    "%8.3lf%8.3lf%8.3lf%8.3lf%8.3lf\n",
                    (double)species_ptr[plant_ptr->sp_idx].cactos_sp_code, 
			        plant_ptr->dbh,
			        plant_ptr->tht,
			        plant_ptr->cr,
			        plant_ptr->expf / n_points );
    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

}


/* MOD003 */
/****************************************************************************/
/*  13. write_cactos_ingrowth_file                                          */
/****************************************************************************/
void __stdcall write_cactos_ingrowth_file(
    unsigned long           *return_code,
    const char              *filename, 
    unsigned long           n_points,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,               /* MOD021   */
    struct COEFFS_RECORD    *coeffs_ptr )           /* MOD021   */
{

	FILE	                *fp;
    struct PLANT_RECORD     *plant_ptr;
    unsigned long           i;
    unsigned long           total_tree_records;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    /* MOD021 */
    /* get the count of the tree records in tha plant list */
    total_tree_records = 0;
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        if( is_tree( &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx] ) )
        {
            total_tree_records++;
        }

    }

    /* export the first line of the file */
	fprintf(    fp, 
/*                 "%-20s%10i\n",  */
                "%-20s%10ld\n",
                filename, 
                total_tree_records );   /* MOD021 */

    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        if( !is_tree( &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx] ) )
        {
            continue;
        }


        fprintf(    fp, 
                    "%8.3lf%8.3lf%8.3lf%8.3lf%8.3lf\n",
                    (double)species_ptr[plant_ptr->sp_idx].cactos_sp_code, 
			        plant_ptr->dbh,
			        plant_ptr->tht,
			        plant_ptr->cr,
			        plant_ptr->expf / n_points );
    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;



}


/* MOD011 */
/****************************************************************************/
/*  14. write_fvs_file                                                      */
/****************************************************************************/
/* this function will iterate through the plots and export the file to an   */
/* FVS file format to be read in by suppose/fvs                             */
/****************************************************************************/
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
    struct COEFFS_RECORD    *coeffs_ptr    )
{

	FILE	                *fp;
    struct PLANT_RECORD     *plant_ptr;
    struct PLOT_RECORD      *plot_ptr;
    unsigned long           i;
    unsigned long           p;
    unsigned long           fvs_cr;
    unsigned long           start_idx;
    unsigned long           end_idx;
    unsigned long           n_plant_recs_on_plot;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    plot_ptr = &plots_ptr[0];
    for( p = 0; p < n_points; p++, plot_ptr++ )
    {
        /* get the plant indecies for the plot */
        get_plant_indecies_for_plot(    return_code,
                                        plot_ptr,
                                        n_plants,
                                        plants_ptr,
                                        &start_idx,
                                        &end_idx,
                                        &n_plant_recs_on_plot );

        /* start writing the plant records to the file */
        plant_ptr = &plants_ptr[start_idx];
        for( i = 0; i < n_plant_recs_on_plot; i++, plant_ptr++ )
        {
            /* you should verify that the species index? for the non-stocked  */
            /* is used instead of checking for "NS"                           */

            /* MOD011 */
            /* MOD012 */
            //if( !is_tree( c_ptr ) && strcmp(plant_ptr->sp_code, "NS") )
            if( !is_tree( &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx] )  && 
                   strcmp( species_ptr[plant_ptr->sp_idx].sp_code, "NS" ) )
                   //strcmp( species_ptr[plant_ptr->sp_idx].sp_code, "NS" ) == 0 )
            {
                continue;
            }

            /* get the prognosis crown ratio */
            fvs_cr = get_fvs_cr( plant_ptr->cr );

            /* get the correct species code from the look up table  */
            /* you should use the variant to find the correct       */
            /* entry from the species table                         */
            /* MOD012  */
            //if( !strcmp(plant_ptr->sp_code, "NS"))
            if( !strcmp( species_ptr[plant_ptr->sp_idx].sp_code, "NS" ) )
            {
                fprintf(    fp,
/*                             "%4i\n", */
                            "%4ld\n",
                            plant_ptr->plot);
            }
            else
            {
	       fprintf(    fp, 
/* 			   "%4i%3i%6.0f1%-3s%4.1f%3s%3.0f%7s%1i\n", */
			   "%4ld%3ld%6.0f1%-3s%4.1f%3s%3.0f%7s%1ld\n",
			   plant_ptr->plot,
			   i+1,
			   plant_ptr->expf,			   
			   species_ptr[plant_ptr->sp_idx].fvs_sp_code,
			   plant_ptr->dbh,
			   " ",
			   plant_ptr->tht,
			   " ",
			   fvs_cr);
            }
        }
    }

    fclose( fp );

    *return_code = CONIFERS_SUCCESS;

}


/****************************************************************************/
/*  15. get_fvs_cr                                                          */
/****************************************************************************/
static int get_fvs_cr( double cr )
{
	int number = 0;

	if( cr == 0.0 ) 
    {
        number = 0;
    }
	if( cr > .01 && cr <= .1 ) 
    {
        number = 1;
    }
	if( cr > .1 && cr <= .2 ) 
    {
        number = 2;
    }
	if( cr > .2 && cr <= .3 ) 
    {
        number = 3;
    }
	if( cr > .3 && cr <= .4 ) 
    {
        number = 4;
    }
	if( cr > .4 && cr <= .5 ) 
    {
        number = 5;
    }
	if( cr > .5 && cr <= .6 ) 
    {
        number = 6;
    }
	if( cr > .6 && cr <= .7 )
    {
        number = 7;
    }
	if( cr > .7 && cr <= .8 )
    {
        number = 8;
    }
	if( cr > .8 && cr <= 1.0 )
    {
        number = 9;
    }

	return number;
}



/*
static int compare_plants_by_plot_and_sp( 
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
        return strcmp( pt1_ptr->sp_code, pt2_ptr->sp_code );
    }

}

*/



/* MOD005 */
/******************************************************************************/
/*  18. combine_plant_lists                                                   */
/******************************************************************************/
/*  Description :   combine_plant_lists                                       */    
/*  Author      :   Jeff D. Hamann                                            */
/*  Date        :   November 13, 1999                                         */
/*  Returns     :   struct PLANT_RECORD *                                     */
/*  Comments    :   function was developed to combine two plant lists         */
/*                  because plant lists may come from different source and    */
/*                  should be combined into a single plant list for each      */
/*                  sample. It is up to the user to free the memory of the    */
/*                  component plant lists that go into the combined plant list*/
/*                  since the new plant list is allocated with the function   */
/*                  and the component plants lists are copied into the        */
/*                  combined plant list.                                      */
/*  Arguments   :                                                             */
/*    unsigned long       *return_code - general return code if mem allocation*/
/*                          was successful and copied the arrays              */
/*    unsigned long       *n_total_plants - pointer the value that hold the   */
/*                          new size of the combined plant lists              */
/*    unsigned long       n_new_plants - number of plants to add to the       */
/*                          combined plant list                               */
/*    struct PLANT_RECORD *original_plants_ptr - pointer to the array of the  */
/*                          original plant list                               */
/*    struct PLANT_RECORD *new_plants_ptr - pointer to the array of the plants*/
/*                          that are to be added to the original plant list   */
/******************************************************************************/
/*  Formula : N/A                                                             */ 
/*  Source  : N/A                                                             */
/*  Coeffs  : N/A                                                             */
/******************************************************************************/
struct PLANT_RECORD *combine_plant_lists( 
    unsigned long       *return_code,
    unsigned long       *n_total_plants,
    unsigned long       n_new_plants,
    struct PLANT_RECORD *original_plants_ptr,
    struct PLANT_RECORD *new_plants_ptr )
{

    struct PLANT_RECORD *temp_plants_ptr;
    
    /* allocate space for both of the plants arrays */
    temp_plants_ptr = 
        (struct PLANT_RECORD *)calloc(
            *n_total_plants + n_new_plants, sizeof( struct PLANT_RECORD ) );

    /* check to make sure you could allocate the new plant list */
    if( temp_plants_ptr == NULL )
    {
        *return_code = CONIFERS_ERROR;
        return NULL;
    }

    /* copy the first plant list into the newly allocated space */
    memcpy( temp_plants_ptr, 
            original_plants_ptr, 
            *n_total_plants * sizeof( struct PLANT_RECORD ) );

    if( temp_plants_ptr == NULL )
    {
        *return_code = CONIFERS_ERROR;
        return NULL;
    }

    /* copy the first plant list into the newly allocated space */
    memcpy( temp_plants_ptr + *n_total_plants, 
            new_plants_ptr, 
            n_new_plants * sizeof( struct PLANT_RECORD ) );

    if( temp_plants_ptr == NULL )
    {
        *return_code = CONIFERS_ERROR;
        return NULL;
    }

    /* now return the reallocated plants array */
    *return_code = CONIFERS_SUCCESS;
    *n_total_plants += n_new_plants;
    return temp_plants_ptr;

}



/*  MOD016 */
/*****************************************************************************/
/*  23. read_systum1_archive                                                 */
/*       returns a SAMPLE_RECORD                                             */
/*       faulty returns:                                                     */
/*          INVALID_FILE_NAME indicates unable to open the file              */
/*          SAMPLE_RECORD_FAILED indicates it couldn't allocate memory       */
/*          CORRUPTED_ARCHIVE    indicates that something was not right with */
/*                                the systum-1 archive file                  */ 
/*****************************************************************************/
void read_systum1_archive(
    unsigned long           *return_code,
    const char              *filename, 
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
	
	unsigned long			*age,
    unsigned long           *n_points,
    struct PLOT_RECORD      **plots_ptr,
    unsigned long           *n_plants,
    struct PLANT_RECORD     **plants_ptr )
{

    struct PLOT_RECORD      *plot_ptr;
    struct PLANT_RECORD     *plant_ptr;

    char                line_buffer[356];
    char                line_buffer2[120];
    unsigned long       i;
    unsigned long       j;
//	char                temp_species[SP_LENGTH];

	/*  temp header variables in the systum-1 archive */
    char                temp_name[20];
    unsigned long       temp_npoints;
    unsigned long       temp_ntrees;
    unsigned long       temp_standage;
    unsigned long       temp_brushage;
    double              temp_ppsi;
    unsigned long       temp_shrub_count;     /* temp count of num shrub records  */

	/*  temp plot variables in the systum-1 archive */
    unsigned long       temp_plot;
    double              temp_sp[6];
    double              temp_cover[6];
    double              temp_ht[6];
    unsigned long       check_count;

	/*  temp plant variables */
    double              plot;
    unsigned long       temp_plant_count;
    double              temp_tree_spcode;
    //char                temp_conif_code[SP_LENGTH];
    double  			temp_dbh;
    double			    temp_tht;
    double			    temp_cr;
    double			    temp_expf;
    double			    temp_pai;

    char                    temp_sp_code[16];
    struct SPECIES_RECORD   *s_ptr;

	FILE	            *fp;

    temp_shrub_count    = 0;
    check_count         = 0;

    if( ( fp = fopen( filename, "rt" ) ) == NULL )
    {
        *return_code = INVALID_FILE_NAME;
        //return NULL;
        return;
    }


	/*----step 1.---------------------------------------------------------------  */    
	/* read the first line of the systum-1 archive file */
    fgets(line_buffer, sizeof( line_buffer ), fp);
    memcpy(temp_name, line_buffer, 19);
    memcpy(line_buffer2, line_buffer + 20 , 59);
    if( sscanf( line_buffer2, "%ld %ld %ld %ld %lf", 
                &temp_npoints,
                &temp_ntrees,
                &temp_standage,
                &temp_brushage,
                &temp_ppsi) < 5)
     {
			fclose( fp );
            *return_code = CORRUPTED_ARCHIVE;
            temp_ntrees  = 0;
            //return NULL;
            return;
     }

    *age			= temp_standage ;
	*n_points		= temp_npoints;

	/*    allocate space for plots                 */
	if( *plots_ptr )
	{
		free( *plots_ptr );
		*plots_ptr = NULL;
	}


    *plots_ptr = (struct PLOT_RECORD*)calloc( *n_points, 
		sizeof( struct PLOT_RECORD ) );


	/*----step 2.-----------------------------------------------------*/
	/*      now find the nubmber of shrub records      */
    plot_ptr = &(*plots_ptr)[0];
    for( i = 0; i < *n_points; i++, plot_ptr++ )
    {
	    fgets( line_buffer, sizeof( line_buffer ), fp );
        
        memcpy( line_buffer2, line_buffer + 54, 106 );
        sscanf(line_buffer, "%ld" , &temp_plot);
//        fgets( line_buffer, sizeof( line_buffer ), fp ); // fix read glitch
        sscanf( line_buffer2, 
		        "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
               " %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
                &temp_sp[0], &temp_cover[0], &temp_ht[0],
                &temp_sp[1], &temp_cover[1], &temp_ht[1],
                &temp_sp[2], &temp_cover[2], &temp_ht[2],
                &temp_sp[3], &temp_cover[3], &temp_ht[3],
                &temp_sp[4], &temp_cover[4], &temp_ht[4],
                &temp_sp[5], &temp_cover[5], &temp_ht[5]);

		/*      this little loop counts the number of nonzero shrubs in the plot */
        for( j = 0; j < 6; j++)
        {
            if( (temp_cover[j] != 0.0) && ( temp_sp[j] != 0.0) ) temp_shrub_count += 1;
		}

        /* initialize the plot values   */
        plot_ptr->plot      = temp_plot;
        plot_ptr->aspect    = 0;
        plot_ptr->elevation = 3000;
        plot_ptr->slope     = 0.0;
		plot_ptr->error		= 0;
        
        /*  MOD019 */
        /*  MOD027 */
		fill_in_whc_and_precip( 
					return_code,
					0.0,
					temp_ppsi,
					&plot_ptr->water_capacity,
					&plot_ptr->mean_annual_precip);
        
        if( *return_code != CONIFERS_SUCCESS )
		{
			plot_ptr->error = *return_code;
		}

    }

    /*----step 3.-------------------------------------------------------------*/
    /*    allocate space for   plant records (trees+shrubs)        */
    *n_plants    = temp_ntrees + temp_shrub_count;
	if( *plants_ptr )
	{
		free( *plants_ptr );
		*plants_ptr = NULL;
	}

    *plants_ptr = (struct PLANT_RECORD *)calloc( 
        (*n_plants), sizeof( struct PLANT_RECORD ) );


/*-----step 4.------------------------------------------------------------*/
/*   now rewind and try it again to fill in the values in the tree record */

    rewind ( fp );

/*-----step 5.------------------------------------------------------------*/
/*    re-read the first line, just to get it out of the way               */

    fgets(line_buffer, sizeof( line_buffer ), fp);

/*----step 6.-------------------------------------------------------------*/
/*   now re-read the plot array and create plant records at top of list   */
    temp_plant_count=0;    


    /* sort the species codes based on sp_code */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_sp_code );

    plant_ptr = &(*plants_ptr)[0];
    for( i = 0; i < *n_points; i++)
    {
	    fgets( line_buffer, sizeof( line_buffer ), fp );

        memcpy( line_buffer2, line_buffer + 54, 104 );
        sscanf( line_buffer, "%ld" , &temp_plot);

        sscanf( line_buffer2, 
		        "%lf %lf %lf %lf %lf %lf %lf %lf %lf" 
                " %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
                &temp_sp[0], &temp_cover[0], &temp_ht[0],
                &temp_sp[1], &temp_cover[1], &temp_ht[1],
                &temp_sp[2], &temp_cover[2], &temp_ht[2],
                &temp_sp[3], &temp_cover[3], &temp_ht[3],
                &temp_sp[4], &temp_cover[4], &temp_ht[4],
                &temp_sp[5], &temp_cover[5], &temp_ht[5]);

        for( j=0; j<6; j++)
        {
            if( ( temp_cover[j] != 0.0 ) && ( temp_sp[j] != 0.0 ) )
            {
	            plant_ptr->plot = temp_plot;
				convert_systum_shrub_species(return_code, 
                                             temp_sp[j],
                                             temp_sp_code );

                /* convert the temporary species code to the internal functional species code */
                s_ptr = get_species_entry_from_code(    n_species,
                                                        species_ptr, 
                                                        temp_sp_code );
                if( !s_ptr )
                {
					*return_code = CONIFERS_ERROR;

					free( *plots_ptr );
					free( *plants_ptr );
					*plots_ptr = NULL;
					*plants_ptr = NULL;

					fclose( fp );
					
					return;
                }

                plant_ptr->sp_idx = s_ptr->idx;

				plant_ptr->pct_cover  = temp_cover[j]*100.0;
				plant_ptr->tht        = temp_ht[j];
                plant_ptr++;
                temp_plant_count++;
            }
        }
    }
    /* next double check to make sure we converted the right number of shrubs */
    if( temp_plant_count != temp_shrub_count )
    {
		free( *plots_ptr );
		free( *plants_ptr );

        *return_code=SYS_ARCHIVE_FAILED;
        return;
    }

	/*----step 7.---------------------------------------------------------------*/        
	/* read rest of the x-array from the systum-1 archive file  */
    for( i = 0; i < temp_ntrees; i++, plant_ptr++ )
    {
    	fgets( line_buffer, sizeof( line_buffer ), fp );

        if( sscanf( line_buffer, 
		        "%lf %lf %lf %lf %lf %lf %lf", 
			      &plot,
			      &temp_tree_spcode,
			      &temp_dbh,
			      &temp_tht,
			      &temp_cr,
			      &temp_expf,
			      &temp_pai ) < 6 )
        {
            continue;/* error trap here */
        }

        plant_ptr->plot     = (int) plot;


        /* is this needed? */
        convert_systum_tree_species(	return_code,
									   temp_tree_spcode,
									   temp_sp_code );
	
	
		/* you need to verify that the species are being converted properly */
		/* convert the temporary species code to the internal functional species code */
		s_ptr = get_species_entry_from_code(    n_species,
							species_ptr, 
							temp_sp_code );
		if( !s_ptr )
		{
		   
		   *return_code = CONIFERS_ERROR;
		   
			free( *plots_ptr );
			free( *plants_ptr );
			*plots_ptr = NULL;
			*plants_ptr = NULL;
		   
		   fclose( fp );
		   return;
		}
		
		plant_ptr->sp_idx = s_ptr->idx;
		
		plant_ptr->dbh      = temp_dbh;
		plant_ptr->tht      = temp_tht;
		plant_ptr->cr       = temp_cr;
		plant_ptr->expf     = temp_expf;

    }
    
    fclose( fp );

    /* then resort the species codes based on idx */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_idx );

	return;

}



/* MOD017 */
/*----------------------------------------------------------------------------*/
/*   24. convert_systum_tree_species                                          */
/*----------------------------------------------------------------------------*/
static void convert_systum_tree_species( 
    unsigned long           *return_code,
    double                  d_code,
    char                    *spcode)
{
//    char    spcode[SP_LENGTH];
    int systum_sp_code;

/*  you must first convert the systum1 code from float to int */

    systum_sp_code = (int) d_code;

    *return_code = CONIFERS_SUCCESS;

    switch (systum_sp_code)
    {
        case 1:
            strcpy(spcode,"PP");
            break;
        case 2:
            strcpy(spcode,"SP");
            break;
        case 3:
            strcpy(spcode,"IC");
            break;
        case 4:
            strcpy(spcode,"DF");
            break;
        case 5:
            strcpy(spcode,"WF");
            break;
        case 6:
            strcpy(spcode,"RF");
            break;
        case 7:
            strcpy(spcode,"LP");
            break;
        case 8:
            strcpy(spcode,"WP");
            break;
        case 9:
            strcpy(spcode,"JP");
            break;
        case 10:
            strcpy(spcode,"OC");
            break;
        case 11:
            strcpy(spcode,"GC");
            break;
        case 12:
            strcpy(spcode,"BO");
            break;
        case 13:
            strcpy(spcode,"TO");
            break;
        case 14:
            strcpy(spcode,"OH");
            break;
        case 15:
            strcpy(spcode,"GF");
            break;
        case 16:
            strcpy(spcode,"WH");
            break;
        case 17:
            strcpy(spcode,"BM");
            break;
        case 18:
            strcpy(spcode,"PM");
            break;
        case 19:
            strcpy(spcode,"LO");
            break;
        case 20:
            strcpy(spcode,"WO");
            break;
        case 21:
            strcpy(spcode,"OC");
            break;
        default:
            *return_code=INVALID_SP_CODE;
            strcpy(spcode,"   ");
    }
}

/* MOD018 */
/*----------------------------------------------------------------------------*/
/*   25. convert_systum_shrub_species                                         */
/*----------------------------------------------------------------------------*/
static void convert_systum_shrub_species( 
    unsigned long           *return_code,
    double                  d_code,
    char                    *spcode)
{

    int systum_shrub_code;

    /*  you must first convert the systum1 code from float to int */
    systum_shrub_code = (int) d_code;

    *return_code = CONIFERS_SUCCESS;

    switch (systum_shrub_code)
    {
        case 1:
            strcpy(spcode,"VM");
            break;
        case 2:
            strcpy(spcode,"BM");
            break;
        case 3:
            strcpy(spcode,"AMSP");
            break;
        case 4:
            strcpy(spcode,"AMSP");
            break;
        case 5:
            strcpy(spcode,"PM");
            break;
        case 6:
            strcpy(spcode,"ARPA");
            break;
        case 7:
            strcpy(spcode,"ARVI");
            break;
        case 8:
            strcpy(spcode,"GC");
            break;
        case 9:
            strcpy(spcode,"CASE3");
            break;
        case 10:
            strcpy(spcode,"CECO2");
            break;
        case 11:
            strcpy(spcode,"CEIN");
            break;
        case 12:
            strcpy(spcode,"CEVE");
            break;
        case 13:
            strcpy(spcode,"CEPR");
            break;
        case 14:
            strcpy(spcode,"CHFO2");
            break;
        case 15:
            strcpy(spcode,"CYSC");
            break;
        case 16:
            strcpy(spcode,"UNKB");
            break;
        case 17:
            strcpy(spcode,"HODI");
            break;
        case 18:
            strcpy(spcode,"TO");
            break;
        case 19:
            strcpy(spcode,"LO");
            break;
        case 20:
            strcpy(spcode,"OO");
            break;
        case 21:
            strcpy(spcode,"CO");
            break;
        case 22:
            strcpy(spcode,"RHDI");
            break;
        case 23:
            strcpy(spcode,"RHMA");
            break;
        case 24:
            strcpy(spcode,"RHOC");
            break;
        case 25:
            strcpy(spcode,"AMSP");
            break;
        case 26:
            strcpy(spcode,"ARPA");
            break;
        case 27:
            strcpy(spcode,"CEVE");
            break;
        case 28:
            strcpy(spcode,"UNKB");
            break;
        case 29:
            strcpy(spcode,"OO");
            break;
        case 30:
            strcpy(spcode,"RHMA");
            break;
        case 31:
            strcpy(spcode,"RISP");
            break;
        case 32:
            strcpy(spcode,"UNKB");
            break;
        case 33:
            strcpy(spcode,"OH");
            break;
        case 34:
            strcpy(spcode,"FORB");
            break;
        case 35:
            strcpy(spcode,"GRAS");
            break;
        case 36:
            strcpy(spcode,"PM");
            break;
        case 37:
            strcpy(spcode,"ARPA");
            break;
        default:
            *return_code=INVALID_SP_CODE;
            strcpy(spcode,"   ");
    }
}



/********************************************************************************/
/* write_organon_inp                                                            */
/********************************************************************************/
void __stdcall write_organon_inp( 
    unsigned long           *return_code,
    const char              *filename, 
    const char              *title, 
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
    struct COEFFS_RECORD    *coeffs_ptr    )
{

    //struct PLOT_RECORD      *plot_ptr;
    struct PLANT_RECORD     *plant_ptr;
//    struct SPECIES_RECORD   *sp_ptr;
//    struct COEFFS_RECORD    *c_ptr;
    char                    even;
    unsigned long           i;
    unsigned long           osp;
    unsigned long           total_tree_records;
    unsigned long           big6_count;
    unsigned long           is_big6;
    unsigned long           other_count;
    unsigned long           num_nest_subsamples;
    FILE                    *fp;

    // Map unknown species and fill in the species group
    //int big6 = sppmap( &tree[0],numoftrees );
    //int other = numoftrees - big6;
    // Calculate the time to breast height
    //int standage = round( StandAge(plot) );
    // Set up the indicator for even aged stands
    // Count the number of points in the tree list
    //int npts = CountPoints(&tree[0],numoftrees);

    /* initialize the variables */
    total_tree_records  = 0;
    big6_count          = 0;
    num_nest_subsamples = 0;

    if( even_age == 0 )
    {
        even = 'F';
    }
    else
    {
        even = 'T';
    }

    /* open the file */
    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_FILE_NAME;
        return;
    }

    /* iterate through the plants and get the count of the big6 records */
    /* write the tree records */
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ ) 
    {
            
        /* if the plant is not a tree, then don't export it */
        if( !is_tree( &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx] ) ) 
        {
            continue;
        }

        if( plant_ptr->expf < 0.0001 )
        {
            continue;
        }

        if( plant_ptr->dbh  < 0.1 )
        {
            continue;
        }

        if( plant_ptr->tht  <= 4.5 )
        {
            continue;
        }

        /* get the organon sp_group */
        osp = get_org_sp_group( return_code,
                                n_species,
                                species_ptr,
                                plant_ptr,
                                &is_big6 );

        total_tree_records++;
        big6_count += is_big6;

    }

    other_count = total_tree_records - big6_count; 

    // Big6 # of tree records in DF,GF,WF,PP,IC,SP
    // Write the header line
    fprintf(    fp,
                "%-12s         %3ld %3ld%5.1lf%5.1lf%1cF%3ld%3ld%4ld%4ld   0",
                title, 
                n_plots, 
                total_tree_records, 
                df_site, 
                pp_site, 
                even, 
                stand_age, 
                bh_age, 
                big6_count, 
                other_count );

    /* 3 = number of nested subsamples      */
    /* 1 = version number                   */
    /* 1 = swo, 2 = nwo, 3 = smc            */
    /* %5.1lf for max sdi's for df, pp, wh  */
    /* in this case the number of nested    */
    /* sub samples is zero                  */
    fprintf(    fp,
                "%1ld%1ld                       0.   0.   0.\n", 
                num_nest_subsamples,
                sim_version );

    /* write the 3 calibration lines */
    for (i = 0; i < 3; i++)
    {
        fprintf(    fp,
                    " 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00\n");
    }

    /* description of nested subsamples */
    fprintf(    fp,
                "     00000  .0  .0  .0  .0  .0    .0    .0    .0    .0    .0\n");

    /* write the tree records */
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ ) 
    {
        /* if the plant is not a tree, then don't export it */
        if( !is_tree( &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx] ) ) 
        {
            continue;
        }

        if( plant_ptr->expf < 0.0001 )
        {
            continue;
        }

        if( plant_ptr->dbh  < 0.1 )
        {
            continue;
        }

        if( plant_ptr->tht  <= 4.5 )
        {
            continue;
        }


        /* get the organon sp_group */
        osp = get_org_sp_group( return_code,
                                n_species,
                                species_ptr,
                                plant_ptr,
                                &is_big6 );

        fprintf(    fp,
                    "%3ld%2ld 0%5.1lf%5.1lf%4.2lf%8.4lf%5.2lf%3ld\n",
                    species_ptr[plant_ptr->sp_idx].organon_sp_code,  
                    osp,
                    plant_ptr->dbh, 
                    plant_ptr->tht, 
                    plant_ptr->cr,
                    plant_ptr->expf,
                    plant_ptr->dbh_growth * 2.5, /* 5-year radial growth */
                    plant_ptr->plot );
    }

    fclose( fp );

}

/********************************************************************************/
/* get_org_sp_group                                                             */
/********************************************************************************/
static unsigned long get_org_sp_group( 
    unsigned long           *return_code,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    struct PLANT_RECORD     *plant_ptr,
    unsigned long           *is_big6 )
{

    
    unsigned long           org_sp_group;
    
    /* set this value first */
    *is_big6         = 0;

    switch( species_ptr[plant_ptr->sp_idx].organon_sp_code )
    {

        case 15: 
        case 19: 
        case 21: 
        case 22: 
        {
            org_sp_group    = 2;
            *is_big6        = 1;
            break;
        }

        case 17: 
        {
            org_sp_group    = 2;
            *is_big6        = 1;
            break;
        }
    
        case 41: 
        case 81: 
        case 242:  
        {
            org_sp_group    = 5;
            *is_big6        = 1;
            break;
        }

        case 103: 
        case 117: 
        case 119: 
        {
            org_sp_group    = 4;
            *is_big6        = 1;
            break;
        }

        case 122: 
        case 116: 
        {
            org_sp_group    = 3;
            *is_big6        = 1;
            break;
        }

        case 202: 
        {
            org_sp_group    = 1;
            *is_big6        = 1;
            break;
        }

        case 263: 
        case 231: 
        {
            org_sp_group    = 6;
            *is_big6        = 1;
            break;
        }

        case 312: 
        case 351: 
        case 352: 
        case 542:  
        {
            org_sp_group = 11;
            break;
        }

        case 361:  
        {
            org_sp_group = 7;
            break;
        }

        case 431: 
        case 981: 
        {
            org_sp_group = 8;
            break;
        }

        case 631:  
        {
            org_sp_group = 9;
            break;
        }

        case 805: 
        case 492: 
        case 760: 
        case 920: 
        case 999:
        {
            org_sp_group = 10;
            break;
        }

        case 815:  
        {
            org_sp_group = 12;
            break;
        }

        case 818:  
        {
            org_sp_group = 13;
            break;
        }
        
        default:
        {
            org_sp_group = 1;
            break;
        }

    }

    return org_sp_group;

}





/* This function was updated to match the current files that are used in the R package */
struct SPECIES_RECORD *read_species_file(
    unsigned long   *return_code,
    const char      *filename, 
    unsigned long   *n_records )
{
    
    #define FIELDS 16 /* MOD003 */
    //#define FIELDS 9 /* MOD003 */

    FILE                    *fp;
    char                    line_buffer[512];
    char                    temp_string[FIELDS][NAME_LENGTH];
    int                     i = 0;
    int                     t = 0;
    char                    seps[]   = ",\t\"";   
    char                    *token;
    struct SPECIES_RECORD   *species_ptr;

//#define SP_LENGTH           6           /*  max length of sp_codes          */   
//#define FVS_SP_LENGTH       3           /*  max length of sp_codes for fvs  */   
//#define NAME_LENGTH         40          /*  max length for long names       */   

	/* buffers for the fields */
	unsigned long		flds = 0;


    species_ptr = NULL;

    *return_code = 0;

    if( ( fp = fopen( filename,"rt" ) ) == NULL )
    {
        *return_code = 1;
        return NULL;
    }

	i = 0;

	/* eat the first (header) line */
    fgets( line_buffer, sizeof( line_buffer ), fp );

    /* count the full lines in the file....                 */
    /* read each product for the rest of the file           */
    /* for now don't do ANY error checking in the line...   */
    while( fgets( line_buffer, sizeof( line_buffer ), fp ) != NULL )
    {

       memset( temp_string, 0, sizeof( temp_string ) );
       
		/* establish string and get first token...  */
        t = 0;
        token = strtok( line_buffer, seps );
        strcpy( temp_string[t], token );

		flds = 1;
        for( t = 1; t < FIELDS; t++ )
        {
	       token = strtok( NULL, seps );
	   
	        if( token != NULL )
	        {
	            strcpy( temp_string[t], token );
			flds++;
		}
        }

		if( flds != FIELDS )
		{
			*return_code = 1;
			fclose( fp );
			return NULL;
		}

		i++;
    }
    
    (*n_records) = i;
    species_ptr = (struct SPECIES_RECORD *)calloc( 
       (*n_records), sizeof( struct SPECIES_RECORD ) );
    
    /* go to the beginning of the file */
    rewind( fp );
  	
	/* eat the first (header) line */
    fgets( line_buffer, sizeof( line_buffer ), fp );

    /* read each product for the rest of the file */
    i = 0;
    while( fgets( line_buffer, sizeof( line_buffer ), fp ) != NULL )
    {
       memset( temp_string, 0, sizeof( temp_string ) );
       
       //memset( &sp1, 0, sizeof( struct SPECIES_RECORD ) );
       
       /* establish string and get first token...  */
        t = 0;
        token = strtok( line_buffer, seps );
        strcpy( temp_string[t], token );

        for( t = 1; t < FIELDS; t++ )
        {
	       token = strtok( NULL, seps );
	   
	        if( token != NULL )
	        {
	            strcpy( temp_string[t], token );
	        }	   
        }

        species_ptr[i].idx			= atoi( temp_string[0] );
        strcpy( species_ptr[i].sp_code, temp_string[1] );
        species_ptr[i].fsp_idx			= atoi( temp_string[2] );
        strcpy( species_ptr[i].common_name,     temp_string[3] );
        species_ptr[i].organon_sp_code		= atoi( temp_string[4] );
        species_ptr[i].cactos_sp_code		= atoi( temp_string[5] );
        strcpy( species_ptr[i].fvs_sp_code,     temp_string[6] );
        species_ptr[i].endemic_mortality	= atof( temp_string[7] );
        species_ptr[i].max_sdi			= atof( temp_string[8] );
        species_ptr[i].browse_damage		= atof( temp_string[9] );
        species_ptr[i].mechanical_damage	= atof( temp_string[10] );
        species_ptr[i].genetic_worth_h		= atof( temp_string[11] );
        species_ptr[i].genetic_worth_d      = atof( temp_string[12] );

		/* added for the swohybrid variant CONIFERS_SWOHYBRID */
		/* climate model coeffs */
        species_ptr[i].min_temp				= atof( temp_string[13] );
        species_ptr[i].max_temp				= atof( temp_string[14] );
        species_ptr[i].opt_temp				= atof( temp_string[15] );


       i++;
    }

    fclose( fp );

    /* now sort the list based on the species index  */
    qsort(  (void*)species_ptr, 
            (size_t)(*n_records), 
            sizeof( struct SPECIES_RECORD ),
	          compare_species_by_idx );

    *return_code = 0;

    return species_ptr;

}



void write_species_file(    
    unsigned long           *return_code,
    const char              *filename, 
    unsigned long           n_records,
    struct SPECIES_RECORD   *species_ptr )
{

	FILE	                *fp;
    struct SPECIES_RECORD   *s;
    unsigned long           i;

    if( ( fp = fopen( filename, "wt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }

    s = &species_ptr[0];
    for( i = 0; i < n_records; i++, s++ )
    {
        fprintf(    fp, 
                    "%ld,\"%s\",%ld,\"%s\",%ld,%ld,\"%s\",%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                    s->idx,
                    s->sp_code,
                    s->fsp_idx,
                    s->common_name,
                    s->organon_sp_code,
                    s->cactos_sp_code,
                    s->fvs_sp_code,
		            s->endemic_mortality,
			        s->max_sdi,
			        s->browse_damage,
			        s->mechanical_damage,
			        s->genetic_worth_h,
                    s->genetic_worth_d,

					/* added for the swohybrid variant CONIFERS_SWOHYBRID */
					s->min_temp,
					s->max_temp,
					s->opt_temp 

					);

    }

    fclose( fp );
    *return_code = CONIFERS_SUCCESS;
}







/****************************************************************************/
/* read_plots_from_file												*/
/* file mimics the read_plots_from_file function							*/
/****************************************************************************/
void read_plants_from_file( 
	unsigned long			*return_code,
    const char				*filename, 
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long			*n_plants,
	struct PLANT_RECORD		**plants_ptr )
{

    char    line_buffer[256];
	unsigned long    i;
	FILE	*fp;
	struct PLANT_RECORD		*plant_ptr = NULL;
    int                     n_args;
    char                    temp_sp_code[16];
    struct SPECIES_RECORD   *s_ptr;


    if( ( fp = fopen( filename, "rt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }


    /* sort the species codes based on sp_code */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_sp_code );

	/* eat the header line */
	fgets( line_buffer, sizeof( line_buffer ), fp );

	i = 0;
	while( fgets( line_buffer, sizeof( line_buffer ), fp ) != NULL )
	{
		i++;
	}

    (*n_plants) = i;
    *plants_ptr = (struct PLANT_RECORD *)calloc( 
        (*n_plants), sizeof( struct PLANT_RECORD ) );

	/* rewind the file pointer to the beginning of the file */
    rewind( fp );

	/* eat the header line */
	fgets( line_buffer, sizeof( line_buffer ), fp );

    plant_ptr = &(*plants_ptr)[0];
    for( i = 0; i < (*n_plants); i++, plant_ptr++ )
    {

      n_args = 0;
      memset( line_buffer, 0, sizeof( line_buffer ) );
      memset( temp_sp_code, 0, sizeof( temp_sp_code ) );

	    fgets( line_buffer, sizeof( line_buffer ), fp );
        
/*
plot sp.code d6 dbh tht cr n.stems expf crown.width
1 DF 1.9 0.2 5.2 0.89 1 100 4.5

*/


        n_args = sscanf( line_buffer, 
            //"%ld %ld %s %lf %lf %lf %lf %lf %lf %ld %lf %lf %lf %ld", 
            "%ld %s %lf %lf %lf %lf %ld %lf %lf", 
            &plant_ptr->plot,           
            //&plant_ptr->plant,          
            temp_sp_code,
            &plant_ptr->d6,             
            //&plant_ptr->d6_area,        
            &plant_ptr->dbh,            
            //&plant_ptr->basal_area,     
            &plant_ptr->tht,            
            &plant_ptr->cr,			    
            &plant_ptr->n_stems,        
            &plant_ptr->expf,           
            &plant_ptr->crown_width );
            //&plant_ptr->crown_area,     
            //&plant_ptr->user_code  );

        /* assign a plant id for accounting */
        /* not required */
        plant_ptr->plant = i;

        /* convert the temporary species code to the internal functional species code */
        //tree_ptr[i].sp_idx = get_sp_idx( temp_sp_code );
        s_ptr = NULL;
        s_ptr = get_species_entry_from_code(    n_species,
                                                species_ptr, 
                                                temp_sp_code );


        /* if the species is invalid for the current species list */
        /* clean up the sample and return null, setting the       */
        /* rc to invalid file                                     */
        if( !s_ptr || n_args < 9 )
        {

            *return_code = CONIFERS_ERROR;

			if( *plants_ptr )
			{
				free( *plants_ptr );
                *plants_ptr = NULL;
			}

            fclose( fp );
			return;

        }

        plant_ptr->sp_idx = s_ptr->idx;

    }


    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

    /* then resort the species codes based on idx */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_idx );



}





/****************************************************************************/
/* read_inp_file												            */
/* file mimics the read_plots_from_file function							*/
/****************************************************************************/
void read_inp_file( 
	unsigned long			*return_code,
    const char				*filename, 
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long			*n_plants,
	struct PLANT_RECORD		**plants_ptr )
{

    char    line_buffer[256];
	unsigned long    i;
	FILE	*fp;
	struct PLANT_RECORD		*plant_ptr = NULL;
    int                     n_args;
    char                    temp_sp_code[16];
    struct SPECIES_RECORD   *s_ptr;


    if( ( fp = fopen( filename, "rt" ) ) == NULL )
    {
        *return_code = INVALID_OPTION;
        return;
    }


    /* sort the species codes based on sp_code */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_sp_code );

    /* eat the four header line for now */
	fgets( line_buffer, sizeof( line_buffer ), fp );
	fgets( line_buffer, sizeof( line_buffer ), fp );
	fgets( line_buffer, sizeof( line_buffer ), fp );
	fgets( line_buffer, sizeof( line_buffer ), fp );
	fgets( line_buffer, sizeof( line_buffer ), fp );

    /* count the number of tree records in the file */
	i = 0;
	while( fgets( line_buffer, sizeof( line_buffer ), fp ) != NULL )
	{
		i++;
	}

    (*n_plants) = i;
    *plants_ptr = (struct PLANT_RECORD *)calloc( 
        (*n_plants), sizeof( struct PLANT_RECORD ) );

	/* rewind the file pointer to the beginning of the file */
    rewind( fp );

	/* eat the header line */
    /* eat the four header line for now */
	fgets( line_buffer, sizeof( line_buffer ), fp );
	fgets( line_buffer, sizeof( line_buffer ), fp );
	fgets( line_buffer, sizeof( line_buffer ), fp );
	fgets( line_buffer, sizeof( line_buffer ), fp );
	fgets( line_buffer, sizeof( line_buffer ), fp );


    plant_ptr = &(*plants_ptr)[0];
    for( i = 0; i < (*n_plants); i++, plant_ptr++ )
    {

      n_args = 0;
      memset( line_buffer, 0, sizeof( line_buffer ) );
      memset( temp_sp_code, 0, sizeof( temp_sp_code ) );

	  fgets( line_buffer, sizeof( line_buffer ), fp );

        /* parse out the fields */
        /* another day of reinventing the wheel... */

        

/*
plot sp.code d6 dbh tht cr n.stems expf crown.width
1 DF 1.9 0.2 5.2 0.89 1 100 4.5

*/


        n_args = sscanf( line_buffer, 
            //"%ld %ld %s %lf %lf %lf %lf %lf %lf %ld %lf %lf %lf %ld", 
            "%ld %s %lf %lf %lf %lf %ld %lf %lf", 
            &plant_ptr->plot,           
            //&plant_ptr->plant,          
            temp_sp_code,
            &plant_ptr->d6,             
            //&plant_ptr->d6_area,        
            &plant_ptr->dbh,            
            //&plant_ptr->basal_area,     
            &plant_ptr->tht,            
            &plant_ptr->cr,			    
            &plant_ptr->n_stems,        
            &plant_ptr->expf,           
            &plant_ptr->crown_width );
            //&plant_ptr->crown_area,     
            //&plant_ptr->user_code  );

        /* convert the temporary species code to the internal functional species code */
        //tree_ptr[i].sp_idx = get_sp_idx( temp_sp_code );
        s_ptr = NULL;
        s_ptr = get_species_entry_from_code(    n_species,
                                                species_ptr, 
                                                temp_sp_code );


        /* if the species is invalid for the current species list */
        /* clean up the sample and return null, setting the       */
        /* rc to invalid file                                     */
        if( !s_ptr || n_args < 9 )
        {

            *return_code = CONIFERS_ERROR;

			if( *plants_ptr )
			{
				free( *plants_ptr );
                *plants_ptr = NULL;
			}

            fclose( fp );
			return;

        }

        plant_ptr->sp_idx = s_ptr->idx;

    }


    fclose( fp );
    *return_code = CONIFERS_SUCCESS;

    /* then resort the species codes based on idx */
    qsort(  (void*)species_ptr, 
            (size_t)(n_species), 
            sizeof( struct SPECIES_RECORD ),
	        compare_species_by_idx );



}





/* this function prints the plant list variables and the associated errors  */
/* it needs to be updated.                                                  */
/* todo: update the function so that it can be more informative             */
//void print_errors_and_warnings(
//    unsigned long           n_plants,
//    struct PLANT_RECORD     *plants_ptr )
//{
//
//
//    unsigned long           i;
//    struct PLANT_RECORD     *p;
//
//    p = &plants_ptr[0];
//    for( i = 0; i < n_plants; i++, p++ )
//    {
//        fprintf(    stdout, 
//                    "plant %ld contains %ld errors\n",
//                    i,
//                    p->errors );
//
//    }
//   
//}



