/****************************************************************************/
/*                                                                          */
/*  sample.c                                                                */
/*  functions used to replicate plots                                       */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                          Revision History                                */
/*                                                                          */
/*  Number  Date          Who     Revision Notes                            */
/****************************************************************************/
/*  MOD000  April 18,2000 MWR     created file and header information       */
/*  MOD001  April 24,2000 MWR     added normal r.deviate and uniform r.d.   */
/*  MOD002  May   11,2000 MWR     tested and modified the uniform r.d. ok   */
/*  MOD003  June  16,2000 MWR     changed n_dups so that you get the        */
/*                                  right num                               */
/*  MOD004  Nov   21,2000 JDH     updated the generate_duplicate_XXX()'s    */
/*  MOD005  Dec   20,2000 JDH     removed static from gaus_dev() since it's */
/*                                  called in other functions and removed   */
/*                                  declaration into conifers.h             */
/****************************************************************************/


/* 	$Id: sample.c 840 2011-11-22 23:47:23Z hamannj $	 */

/* #ifndef lint */
/* static char vcid[] = "$Id: sample.c 840 2011-11-22 23:47:23Z hamannj $"; */
/* #endif /\* lint *\/ */


//#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "conifers.h"


/****************************************************************************/
/*  calc_replication_factor                                                 */
/****************************************************************************/
/*  Description :   performs basic error checking for the samnple record    */   
/*  Author      :   Martin W. Ritchie                                       */
/*  Date        :   April 18, 2000                                          */
/*  Returns     :   void                                                    */
/*  Comments    :   this is a function to calculate how many dup. plots     */
/*                  to make                                                 */
/*  Arguments   :   *return_code                                            */
/*                  n_points            number of points read in the sample */
/*                  n_plants            number of plants read in the sample */
/*  return int                                                              */
/****************************************************************************/
unsigned long calc_replication_factor(
    unsigned long           *return_code,
    unsigned long           n_points,
    unsigned long           n_plants,
    unsigned long           max_sample_size )
{
    double                  temp_replication;
    unsigned long           replication_factor;

    temp_replication = 0.0;
    replication_factor = 0;
    if(max_sample_size ==0)
    {
        replication_factor=0;
        return replication_factor;
    }

    temp_replication   = ( (double)max_sample_size / (double)n_plants );
    
    replication_factor = (unsigned long)temp_replication;

    return replication_factor;

}



/* MOD004   */
/* MOD007   */
/****************************************************************************/
/* generate_duplicate_plots                                                 */
/****************************************************************************/
/*  Description :   generates n duplicates for each plot                    */
/*  Author      :   Jeff D. Hamann/Martin W. Ritchie                        */
/*  Date        :   May 23, 2000                                            */
/*  Returns     :   void                                                    */
/*  Comments    :   this is a utility function that copies the              */
/*                  data from the original plot list and generates          */
/*                  replicates from the original data                       */
/*  Arguments   :                                                           */
/*  unsigned long *return_code  -   return code for calling function        */
/*                  to check                                                */
/****************************************************************************/
/*  Formula : NONE                                                          */
/*  Source  : NONE                                                          */
/*  Coeffs  : NONE                                                          */
/****************************************************************************/
struct PLOT_RECORD  *generate_duplicate_plots(
    unsigned long       *return_code, 
    unsigned long       *n_new_plots,
    unsigned long       n_dups,
    unsigned long       n_orig_plots,
    struct PLOT_RECORD  *orig_plots_ptr )
{

    unsigned long       i;
    unsigned long       p;
    unsigned long       max_plot_num;

    struct PLOT_RECORD  *temp_plots_ptr;
    struct PLOT_RECORD  *copy_plot_buffer;  /* buffer for copying */    
    struct PLOT_RECORD  *max_plot_ptr;

    /* allocate the new plots */
    temp_plots_ptr = (struct PLOT_RECORD* )calloc( 
        n_orig_plots + ( n_orig_plots * n_dups ), 
        sizeof( struct PLOT_RECORD ) );
    
    /* create a temp buffer for copying the plots */
    copy_plot_buffer = (struct PLOT_RECORD* )calloc( 
        n_orig_plots, sizeof( struct PLOT_RECORD ) );

    /* copy the original data into the buffer */
    memcpy( copy_plot_buffer, 
            orig_plots_ptr, 
            n_orig_plots * sizeof( struct PLOT_RECORD ) );

    
    /* iterate throught the plots and copy the array    */
    /* into the next block                              */
    for( i = 0; i < n_dups + 1; i++ )
    {
        /* set up the buffer for copying but only do    */
        /* this if you're past the first "duplicate"    */
        max_plot_num = 0;
        max_plot_ptr = &copy_plot_buffer[0];
        for( p = 0; p < n_orig_plots; p++, max_plot_ptr++ )
        {
            if( (unsigned long)max_plot_ptr->plot > max_plot_num )
            {
                max_plot_num = max_plot_ptr->plot;
            }
        }
        
        /* now, apply the max to the buffer that will get copied in */
        if( i > 0 ) 
        {
            for( p = 0; p < n_orig_plots; p++ )
            {
                copy_plot_buffer[p].plot = orig_plots_ptr[p].plot + max_plot_num;
            }
        }

        /* copy the plot data into the correct place    */
        /* in the new plot array                        */
        memcpy( temp_plots_ptr + (n_orig_plots * i), 
                copy_plot_buffer, 
                n_orig_plots * sizeof( struct PLOT_RECORD ) );

    }


    /* free up the temp buffer */
    free( copy_plot_buffer );

    /* update the plots */    
    *n_new_plots    = n_orig_plots * ( n_dups + 1 );
    *return_code    = CONIFERS_SUCCESS;

    return temp_plots_ptr;
}


/* MOD004   */
/* MOD008   */
/****************************************************************************/
/* generate_duplicate_plants                                                */
/****************************************************************************/
/*  Description :   generates n duplicates for each plot                    */
/*  Author      :   Jeff D. Hamann/Martin W. Ritchie                        */
/*  Date        :   May 23, 2000                                            */
/*  Returns     :   void                                                    */
/*  Comments    :   this is a utility function that copies the data         */
/*                  from the what                                           */
/*  Arguments   :                                                           */
/*  unsigned long *return_code  -   return code for calling function        */
/*  unsigned long       n_dups  -   number of duplicate records to be made  */
/*  unsigned long       n_plants -  number of plants in the original list   */
/*  unsigned long       n_total_plots number of plots after duplication     */
/*  struct PLANT_RECORD *plants_ptr  - plants pointer before duplication    */
/*  struct PLOT_RECORD  *plots_ptr   - plots ptr after duplication          */
/****************************************************************************/
/*  Formula : NONE                                                          */
/*  Source  : NONE                                                          */
/*  Coeffs  : NONE                                                          */
/****************************************************************************/
struct PLANT_RECORD  *generate_duplicate_plants(
    unsigned long       *return_code, 
    unsigned long       *n_new_plants,
    unsigned long       n_dups,
    unsigned long       n_orig_plants,
    struct PLANT_RECORD  *orig_plants_ptr )
{

    unsigned long       i;
    unsigned long       p;
    unsigned long       max_plot_num;

    struct PLANT_RECORD  *temp_plants_ptr;
    struct PLANT_RECORD  *copy_plant_buffer;  /* buffer for copying */    
    struct PLANT_RECORD  *max_plant_ptr;

    /* allocate the new plots */
    temp_plants_ptr = (struct PLANT_RECORD* )calloc( 
        n_orig_plants + ( n_orig_plants * n_dups ), 
        sizeof( struct PLANT_RECORD ) );
    
    /* create a temp buffer for copying the plots */
    copy_plant_buffer = (struct PLANT_RECORD* )calloc( 
        n_orig_plants, sizeof( struct PLANT_RECORD ) );

    /* copy the original data into the buffer */
    memcpy( copy_plant_buffer, 
            orig_plants_ptr, 
            n_orig_plants * sizeof( struct PLANT_RECORD ) );

    
    /* iterate throught the plots and copy the array    */
    /* into the next block                              */
    for( i = 0; i < n_dups + 1; i++ )
    {
        /* set up the buffer for copying but only do    */
        /* this if you're past the first "duplicate"    */
        max_plot_num = 0;
        max_plant_ptr = &copy_plant_buffer[0];
        for( p = 0; p < n_orig_plants; p++, max_plant_ptr++ )
        {
            if( (unsigned long)max_plant_ptr->plot > max_plot_num )
            {
                max_plot_num = max_plant_ptr->plot;
            }
        }
        
        /* now, apply the max to the buffer that will get copied in */
        if( i > 0 ) 
        {
            for( p = 0; p < n_orig_plants; p++ )
            {
                copy_plant_buffer[p].plot = orig_plants_ptr[p].plot + max_plot_num;
            }
        }

        /* copy the plot data into the correct place    */
        /* in the new plot array                        */
        memcpy( temp_plants_ptr + (n_orig_plants * i), 
                copy_plant_buffer, 
                n_orig_plants * sizeof( struct PLANT_RECORD ) );

    }


    /* free up the temp buffer */
    free( copy_plant_buffer );

    /* update the plots */    
    *n_new_plants   = n_orig_plants * ( n_dups + 1 );
    *return_code    = CONIFERS_SUCCESS;

    return temp_plants_ptr;
}



/* MOD005   */
/****************************************************************************/
/* gauss_dev                                                                */
/****************************************************************************/
/* this function returns a gaussian deviant                                 */
/* from NRC                                                                 */
/****************************************************************************/
//static float gauss_dev()
float gauss_dev()
{
    double  fac = 0.0;
    double  v1  = 0.0;
    double  v2  = 0.0;
    double  r   = 0.5;

    for(;;)
    {
       //v1 = 2.0f * ( (float)rand() / 32768.0f ) - 1.0f;
       // v2 = 2.0f * ( (float)rand() / 32768.0f ) - 1.0f;
       
       v1 = 2.0f * ( (double)rand() / (double)RAND_MAX ) - 1.0f;
       v2 = 2.0f * ( (double)rand() / (double)RAND_MAX ) - 1.0f;

       r  = v1 * v1 + v2 * v2;

        if(  r > 0.0f && r < 1.0f )
        {
            break;
        }
    }

    fac = sqrt( -2.0f * log( r ) / r );

    return (float)( v2 * fac );
}



/****************************************************************************/
/* uniform_0_1                                                              */
/****************************************************************************/
/* this function generates a uniform distribution number between 0 and 1    */
/****************************************************************************/
float uniform_0_1()
{
    float i1;

//    i1 =  (float)rand() / (float)(RAND_MAX + 1.0);
    i1 =  (float)( (double)rand() / (double)(RAND_MAX + 1.0) );

    return i1;
}


/****************************************************************************/
/* fill in missing tree expfs                                               */
/****************************************************************************/
/* this function generates a uniform distribution number between 0 and 1    */
/****************************************************************************/
void fill_in_missing_tree_expf(
    unsigned long       *return_code, 
    double              fixed_plot_radius,
    double              min_dbh,
    double              baf,
    unsigned long       n_plants,
    struct PLANT_RECORD *plants_ptr )
{

    unsigned long       i;
    struct PLANT_RECORD *plant_ptr;

    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        if( plant_ptr->dbh > min_dbh )
        {
            plant_ptr->expf = baf / (plant_ptr->dbh * plant_ptr->dbh * FC_I);
        }
        else
        {
            plant_ptr->expf = 43560.0 / 
                ( fixed_plot_radius * fixed_plot_radius * MY_PI );
        }            

        /* multiply the expansion factor by the number of   */
        /* stems that this plant record represents          */
        plant_ptr->expf *= plant_ptr->n_stems;
    }

    *return_code = CONIFERS_SUCCESS;
    
}


