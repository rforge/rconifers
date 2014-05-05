
/* 	$Id: coeffs_mgt.c 861 2012-02-22 21:25:32Z hamannj $	 */

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"

/* local function declarations */
static int compare_coeffs_by_fsp( 
    const void *ptr1, 
    const void *ptr2 );


/* this is an interface function that needs to be exposed at least      */
/* when the library is being initialized                                */
/* todo: step #5 - add init_coeffs function for each new variant here   */
struct COEFFS_RECORD *con_init_coeffs( 
	unsigned long				version,
	unsigned long				*n_coeffs,
    double                      *coeffs_version,
    double                      *model_version )
{

    struct COEFFS_RECORD *c;

    *n_coeffs = 0;

   switch( version )
   {
      case CONIFERS_SWO:
	    c = (struct COEFFS_RECORD *)con_swo_init_coeffs( n_coeffs, coeffs_version, model_version );
	    break;

      case CONIFERS_SMC:
	    c = (struct COEFFS_RECORD *)con_smc_init_coeffs( n_coeffs, coeffs_version, model_version );
	    break;

      case CONIFERS_SWOHYBRID:
	    c = (struct COEFFS_RECORD *)con_swo_hybrid_init_coeffs( n_coeffs, coeffs_version, model_version );
	    break;

      case CONIFERS_CIPS:
	    //c = (struct COEFFS_RECORD *)con_cips_init_coeffs( n_coeffs, coeffs_version, model_version );
	    break;

      default:
	    c = (struct COEFFS_RECORD *)con_swo_init_coeffs( n_coeffs, coeffs_version, model_version );
	    break;
   }

    /* now sort the coeffs for the lookup function */
    qsort(  (void*)c, 
          (size_t)(*n_coeffs), 
          sizeof( struct COEFFS_RECORD ),
          compare_coeffs_by_fsp );
    
    return c;

}



/********************************************************************************/
/* qsort driver function                                                        */
/********************************************************************************/
static int compare_coeffs_by_fsp( 
    const void *ptr1, 
    const void *ptr2 )
{
    struct COEFFS_RECORD   *sp1_ptr;
    struct COEFFS_RECORD   *sp2_ptr;

    sp1_ptr = (struct COEFFS_RECORD*)ptr1;
    sp2_ptr = (struct COEFFS_RECORD*)ptr2;

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

/********************************************************************************/
/* is_tree                                                                      */
/********************************************************************************/
unsigned long __stdcall is_tree( 
    struct COEFFS_RECORD *c_ptr )
{

    unsigned long ret_value;

	if( c_ptr->type == CONIFER || c_ptr->type == HARDWOOD )
    {
        ret_value =  1;
    }
    else
    {
        ret_value =  0;
    }

    return ret_value;

}

/********************************************************************************/
/* is_shrub                                                                     */
/********************************************************************************/
unsigned long __stdcall is_shrub( 
    struct COEFFS_RECORD *c_ptr )
{

    unsigned long ret_value;
    if( c_ptr->type == SHRUB )
    {
        ret_value =  1;
    }
    else
    {
        ret_value =  0;
    }

    return ret_value;

}


/********************************************************************************/
/* is_shrub                                                                     */
/********************************************************************************/
unsigned long __stdcall is_forb( 
    struct COEFFS_RECORD *c_ptr )
{

    unsigned long ret_value;
    if( c_ptr->type == FORB )
    {
        ret_value =  1;
    }
    else
    {
        ret_value =  0;
    }

    return ret_value;

}


/********************************************************************************/
/* is_non_stocked                                                               */
/********************************************************************************/

/* is this really needed */
unsigned long __stdcall is_non_stocked( 
    struct COEFFS_RECORD *c_ptr )
{
    unsigned long ret_value;

    if( c_ptr->type == NON_STOCKED )
    {
        ret_value =  1;
    }
    else
    {
        ret_value =  0;
    }

    return ret_value;


}


