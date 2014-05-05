/********************************************************************************/
/*                                                                              */
/*  plot.c                                                                      */
/*  functions used to predict values for the CONIFERS growth model              */
/*                                                                              */
/********************************************************************************/

/********************************************************************************/
/*                          Revision History                                    */
/*                                                                              */
/*  Number  Date        Who     Revision Notes                                  */
/********************************************************************************/
/*  MOD000  Sept12,1999 JDH     Created file and added header information       */
/*  MOD001  Sept14,1999 JDH     Removed #include "model.h" and incorped it      */
/*                              into "usfs_r5.h"                                */
/*  MOD002  Sept20,1999 JDH     fixed warning that all control points don't     */
/*                              return value, now returns pointer that was      */
/*                              allocated                                       */
/*  MOD003  Sept21,1999 JDH     moved  import_plot_array_from_file() and        */
/*                              write_plot_file() from plot.c to file_io.c      */
/*  MOD004  Nov 12,1999 JDH     fixed bug in build_plot_array_from_plants       */
/*                              to increment through the plants_ptr             */
/*  MOD005  Nov 12,1999 JDH     added convert_dd_2_dms() and convert_dms_2_dd   */
/*  MOD006  Dec 06,1999 JDH     added copy_point_data()                         */
/*                                                                              */
/********************************************************************************/

/* 	$Id: plot.c 854 2012-02-06 23:37:24Z hamannj $	 */

/* #ifndef lint */
/* static char vcid[] = "$Id: plot.c 854 2012-02-06 23:37:24Z hamannj $"; */
/* #endif /\* lint *\/ */


//#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"


/* local functions */
static int compare_long( 
    const void *ptr1, 
    const void *ptr2 );

static int compare_plots( 
    const void *ptr1, 
    const void *ptr2 );


/****************************************************************************/
/* implimentation of functions                                              */
/****************************************************************************/

/* this function examines the plants data and creates, in memory, an		*/
/* array of PLOT_RECORDS													*/
struct PLOT_RECORD *build_plot_array_from_plants( 
    unsigned long       *return_code,
    unsigned long       n_records, 
    struct PLANT_RECORD *plants_ptr,
    unsigned long       *n_points )
{

	  unsigned long    *point_array;
	  unsigned long    point_found;
	 unsigned  long    pt;
	 unsigned long    i;
	 unsigned long    j;
	 unsigned long    min_point;
	 unsigned long    max_point;
    struct  PLOT_RECORD     *plot_ptr;
    struct  PLANT_RECORD    *plant_ptr;

    *n_points = 0;
    point_array = (unsigned long *)calloc( n_records, sizeof( unsigned long ) );

    if( point_array == NULL )
    {
        /* couldn't allocate the point_array */
        *return_code    = CONIFERS_ERROR;
        *n_points       = 0;
        return NULL;
    }

    /* find the min and max point values in the array   */
    /* MOD004   */
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < (unsigned long)n_records; i++, plant_ptr++ )
    {
        point_array[i] = plant_ptr->plot;
    }

    /* sort the array from min to max and get the extreams */
    qsort(  (void*)point_array, 
            (size_t)n_records, 
            sizeof( unsigned long ),
		    compare_long );

    min_point = point_array[0];
    max_point = point_array[n_records-1];
    
    /* deallocate the point array, since you don't need it anymore */
    free( point_array );

    /* count the number of unique plots in the tree list            */
    /* from the smallest found to the largest, count the number of  */
    /* points inbetween the two                                     */
    for( i = min_point; i <= max_point; i++ )
    {
        
        point_found = 0;

        /* go though the trees and see if the "next point is    */
        /* in the array                                         */
        /* find the min and max point values in the array       */
        plant_ptr = &plants_ptr[0];
        for( j = 0; j < (unsigned long)n_records; j++, plant_ptr++ )
        {
            /* if you found a match */
            if( plant_ptr->plot == i && !point_found )
            {
                point_found = 1;
                break;
            }
        }

        if( point_found )
        {
            (*n_points)++;
        }

    }
    
    /* generate the plot array now                  */
    /* you'll have to go through the array again    */
    /* cause there's no ansi _msize()               */
    plot_ptr = (struct PLOT_RECORD*)calloc( 
        (*n_points), sizeof( struct PLOT_RECORD ) );

    if( plot_ptr == NULL )
    {
        /* couldn't allocate the point_array */
        *return_code    = CONIFERS_ERROR;
        *n_points       = 0;
        return NULL;
    }

    
    pt = 0;
    for( i = min_point; i <= max_point; i++ )
    {
        
        point_found = 0;

        /* go though the trees and see if the "next point is    */
        /* in the array                                         */
        /* find the min and max point values in the array       */
        plant_ptr = &plants_ptr[0];
        for( j = 0; j < ( unsigned long)n_records; j++, plant_ptr++ )
        {
            /* if you found a match */
            if( plant_ptr->plot == i && !point_found )
            {
                point_found = 1;
                break;
            }
        }

        if( point_found )
        {
            plot_ptr[pt].plot = i;
            pt++;
        }

    }


    /* MOD002 */
    *return_code = CONIFERS_SUCCESS;
    return plot_ptr;
}


/****************************************************************************/
/* sorting functions for species records                                    */
/****************************************************************************/
static int compare_long( 
    const void *ptr1, 
    const void *ptr2 )
{
	 unsigned     long    *sp1_ptr;
	 unsigned     long    *sp2_ptr;

    sp1_ptr = (	 unsigned long*)ptr1;
    sp2_ptr = (	 unsigned long*)ptr2;

    if( *sp1_ptr < *sp2_ptr )
    {
        return -1;
    }
    if( *sp1_ptr > *sp2_ptr )
    {
        return 1;
    }
    else
    {
        return 0;
    }

}


static int compare_plots( 
    const void *ptr1, 
    const void *ptr2 )
{
    struct PLOT_RECORD   *pt1_ptr;
    struct PLOT_RECORD   *pt2_ptr;

    pt1_ptr = (struct PLOT_RECORD*)ptr1;
    pt2_ptr = (struct PLOT_RECORD*)ptr2;

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

struct PLOT_RECORD *get_plot( 
    long                plot,
    unsigned long       n_records,
    struct PLOT_RECORD  *plots_ptr )
{

    struct PLOT_RECORD   *entry;
	struct PLOT_RECORD   key;
	memset( &key, 0, sizeof( struct PLOT_RECORD ) );
    key.plot = plot;

    /* find the matching plot record */
	entry = (struct PLOT_RECORD*)bsearch( 
        &key, 
		(const void*)plots_ptr, 
		(size_t)n_records,
		sizeof( struct PLOT_RECORD ),
		compare_plots );

    return entry;

}


/* this function fills in the starting and ending index values  */
/* with the starting and ending index values that correspond    */
/* to the plot number that is passed into the function          */
/* it's important to note that the plant array needs to be      */
/* sorted by plot before this function is called                */
void get_plant_indecies_for_plot(    
    unsigned long       *return_code,
    struct PLOT_RECORD  *plot_ptr,
    unsigned long       n_plants,
    struct PLANT_RECORD *plants_ptr,
    unsigned long       *start_idx,
    unsigned long       *end_idx,
    unsigned long       *n_plant_records_on_plot )
{

    unsigned long       i;
    int                 first_record_found;
    int                 last_record_found;
    struct PLANT_RECORD *plant_ptr;


    *start_idx                  = 0;
    *end_idx                    = 0;
    *n_plant_records_on_plot    = 0;

    first_record_found  = 0;
    last_record_found   = 0;

    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        
        /* if you've got a plant on the plot    */
        /* add on to the number of plants on    */
        /* the plot                             */
        if( plant_ptr->plot == plot_ptr->plot )
        {
            if( first_record_found == 0 )
            {
                *start_idx          = i;
                first_record_found  = 1;
            }

            (*n_plant_records_on_plot)++;
        }
    
    }    

    /* since end_idx is an unsigned long it is possible that the    */
    /* value can overflow and this should be fixed!                 */
    if( *start_idx + (*n_plant_records_on_plot) > 0 )
    {
        *end_idx = *start_idx + (*n_plant_records_on_plot) - 1;
        *return_code = CONIFERS_SUCCESS;
    }
    else
    {
        *return_code = CONIFERS_ERROR;
    }

}


/* MOD005 */
//void convert_dd_2_dms( 
//    unsigned long   *return_code,
//    double          dec_deg,
//    double          *dd,
//    double          *mm,
//    double          *ss )
//{

//    *dd = 0.0;
//    *mm  =0.0;
//    *ss = 0.0;

    /* get the whole part of the dec_deg for the hours  */
    /* take the remained and pass that into modf to get */
    /* the minutes                                      */
//    *mm = modf( dec_deg, dd );
//    *mm *= 60.0;

    /* convert if negative */
//    if( *mm < 0 )
//    {
//        *mm *= -1.0f;
//    }

    /* do the same for the seconds */
//    *ss = modf( *mm, mm );
//    *ss *= 60.0;

//    *return_code = CONIFERS_SUCCESS;

//}

/* MOD005 */
//void convert_dms_2_dd( 
//    unsigned long   *return_code,
//    double          dd,
//    double          mm,
//    double          ss,
//    double          *dec_deg )
//{
//
//    *dec_deg = 0.0;
//
//    *dec_deg = dd + ( ( mm + ss / 60.0 ) / 60.0 );
//
//    *return_code = CONIFERS_SUCCESS;
//
//}

/* MOD006 */
/********************************************************************************/
/* copy_point_data                                                              */
/********************************************************************************/
/*  Description :   copies the src plot to the rest of the plots in the         */
/*                  plots_ptr   array                                           */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   December 06, 1999                                           */
/*  Returns     :   void                                                        */
/*  Comments    :   this is a utility function that copies the data from the    */
/*                  source point into the rest of the points in the plots_ptr   */
/*                  array.                                                      */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  src_plot                    -   source plot that will be copied to the rest */
/*                                  of the plots in the plots_ptr array         */
/*  n_points                    -   number of elements in the plots_ptr array   */
/*  struct PLOT_RECORD  *plots_ptr  -   pointer to an array of the plots        */
/********************************************************************************/
/*  Formula : NONE                                                              */
/*  Source  : NONE                                                              */
/*  Coeffs  : NONE                                                              */
/********************************************************************************/
void    copy_point_data( 
    unsigned long           *return_code,
    unsigned long           n_points,
    struct PLOT_RECORD      *src_plot_ptr,
    struct PLOT_RECORD      *dest_plots_ptr )
{
    
    /* local variables */    
    struct PLOT_RECORD  *plot_ptr;
    unsigned long       i;
    

    /* just copy the src_plot_ptr into the dest_plots_ptr elements */
    plot_ptr   = &dest_plots_ptr[0];
    for( i = 0; i < n_points; i++, plot_ptr++ )
    {
        memcpy( plot_ptr, src_plot_ptr, sizeof( struct PLOT_RECORD ) );
        plot_ptr->plot = i + 1;
    }


    *return_code = CONIFERS_SUCCESS;

}



/* this function will iterate through the tree list */
/* and reduce the expf's and percent cover values   */
/* for the species code passed into the function    */
void reduce_pct_cover( 
    unsigned long           *return_code,
    double                  target_pct,
    double                  current_pct,
    unsigned long           target_sp, 
    unsigned long           n_plants, 
    unsigned long           n_plots,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plots_ptr)
//    unsigned long           n_species,
//    struct SPECIES_RECORD   *species_ptr )
{

    unsigned long       i;
    unsigned long       j;
    double              expf_redux;
    struct PLANT_RECORD *plant_ptr;
    double              current_cov;
    double              ending_cov;
    struct PLOT_RECORD *plot_ptr;
    double              plot_cover;
    /* don't bother if the thinning will add    */
    /* to the expansion factor                  */
    if( target_pct < 0.0 || target_pct > 100.0 )
    {
        return;
    }

    expf_redux = 0.0;
    
    plant_ptr = &plants_ptr[0];
    for( i = 0; i < n_plants; i++, plant_ptr++ )
    {
        /* you've found a entry that matches the species    */
        /* so reduce the expf and recalc the pct_cover      */
        /* pct_cover value for the record                   */
        if( plant_ptr->sp_idx == target_sp )
        {
            plot_ptr   = &plots_ptr[0];
            plot_cover = 0.0;

            for( j = 0; j < n_plots; j++, plot_ptr++)
            {
                if(plot_ptr->plot == plant_ptr->plot)
                {
                    plot_cover = plot_ptr->shrub_pct_cover;
                }
            }
            current_cov = 100*(MY_PI*plant_ptr->expf*plant_ptr->crown_width*plant_ptr->crown_width/4.0)/SQ_FT_PER_ACRE;
            ending_cov  = 100.0-(double)target_pct;

            if( ending_cov < current_pct && ending_cov >= 0.0 && current_pct >= 0.0)
            {
                plant_ptr->expf *= (ending_cov / current_pct);
            }
        }
    }
}

