
/********************************************************************************/
/*                                                                              */
/*  mortality.c                                                                 */
/*  functions used to do mortality stuff for the CONIFERS growth model          */
/*                                                                              */
/********************************************************************************/

/* 	$Id: mortality.c 930 2014-01-29 21:48:28Z mritchie $	 */

/*------------------------------------------------------------------------------*/
/*  functions index                                                             */
/*    1. calc_sdi_mortality                                                     */
/*    2. calc_hann_wang_xo                                                      */
/*                                                                              */
/*------------------------------------------------------------------------------*/


#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "conifers.h"


/********************************************************************************/
/*                  calc_sdi_mortality                                          */
/********************************************************************************/
/*  Description :   calculate the background mortality                          */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   January 8, 2000                                             */
/*  Returns     :   void                                                        */
/*  Comments    :   This should be called only if user indicates it is ON       */
/*  Arguments   :   *return code                                                */
/*   return void                                                                */
/*   unsigned	long *return_code  - return code for calling function           */
/*   double		qmd -- the stand quadratic mean diameter                        */
/*	 double		sdimax -- the stand wieghted avg max sdi value                  */
/*   double		x0 -- the hann and wang initial x0 for trajectory               */
/*   double		bhtpa -- tpa of trees > 4.5 feet tall                           */
/*   double		*mortality_prpoportion -- proportion to reduce exp              */
/********************************************************************************/
/*  Formula : N/A                                                               */ 
/*  Source  : Hann and Wang 1990 Research Bulletin 67                           */
/*  Coeffs  : imbedded from Hann and Wang                                       */
/********************************************************************************/
void calc_sdi_mortality(
   unsigned long   *return_code,
   double          qmd,
   double          sdi,
   double          sdimax,
   double          x0,
   double          bhtpa,
   double          *mortality_proportion)
{


   /* new variables for rewrite */
   double myi;
   double a1;
   double a2;
   double a3;
   double yi;
   double xi;
   double y0;
   //double rd;
   int    i;
   double qq;			/* binary search interval */
   double tpa_iter;    /* target tpa from binary search */

   *mortality_proportion = 0.0;

   /* trap zero sdimax I am not sure what else to do with this. */
   /* this should never happen because we only call for trees   */
   if( sdimax <= 0.0001 )
   {
      *return_code = CONIFERS_SUCCESS;
      return;
   }

   /* now calculate if the sdi mortality is needed   */
   /* if not, we are out of here                     */
   if ( bhtpa < 0.0 || x0 <= 0.0 )
   {
      *return_code=CONIFERS_SUCCESS;
      *mortality_proportion = 0.0;
      return;
   }


   /*  the default mortality rate is zero: */
   if( ( sdi/sdimax ) < CRITICAL_RD )
   {
      *return_code = CONIFERS_SUCCESS;
      return;
   }

   /*  from Hann and Wang Paper */
   a3 = 1.47343f;
   a2 = -1.0f / REINEKE_B1;
   a1 = log(10.0f) - a2 * log(sdimax);
   xi = log(bhtpa);
   //rd = sdi / sdimax;

   /* Hann and Wang equation 3*/
   myi = a1 + a2 * log(bhtpa);

   /* calculate y0 as a function of x0 and sdimax H&W page 10 */
   y0 = log(10.0f) - a2 * (log( CRITICAL_RD * sdimax) ) + a2 * x0;

   /* caluculate yi H&W page 8, equation 4*/
   yi = myi - (a1 + a2 * x0 - y0) * exp( - a3 * (x0 - xi) );


   /* the .001 correction here is to handle real close comparisons */
   /* or the very first iteration. Basically anything real close   */
   /* to the Hann and Wang trajectory                              */
   if( qmd  <= exp (yi) + .001f )
   {
      *mortality_proportion = 0.0;    /* then don't do any mortality  */
      *return_code = CONIFERS_SUCCESS; /* this is H&W step 3          */
      return;
   }
   else	/*  then kill trees   */
   {
      tpa_iter = bhtpa * 0.5;
      qq = tpa_iter;		

      /*  set up a loop to find the solution ?    */
      /* is this the interval bisection method?   */
      for( i = 0; i < 500000; i++ )
      {
         qq *= 0.5;
         yi = (a1 + a2 * log(tpa_iter)) - (a1 + a2 * x0 - y0) 
	        * exp( - a3 * (x0 - log(tpa_iter)) );

         if( yi > log( qmd ) )
         {
	        tpa_iter += qq;
         }
         else
         {
	        tpa_iter -= qq;
         }

         /* todo: make sure this can pass mutli-platform precision checks   */ 
         //if( fabs( exp( yi ) - qmd ) < 0.0001 )
         //if( fabs( exp( yi ) - qmd ) < 0.1 )
         //if( fabs( exp( yi ) - qmd ) < 0.0000001 )
         //2.220446e-16

         //if( fabs( exp( yi ) - qmd ) < 2.220446e-16 )  /* double precision */
         if( fabs( exp( yi ) - qmd ) < 1.192093e-07 )   /* single precision */
         {
	        *mortality_proportion = (double)(1.0 - tpa_iter / bhtpa);
	        *return_code = CONIFERS_SUCCESS;
	        return;
         }
      }
    
      /*   if you get here we have a problem   */
      *return_code = CONIFERS_ERROR;
      *mortality_proportion = 1.0f;
      return;
   } /* end of else-kill trees code */

}




/*  MOD020   */
/*  MOD023   */
/********************************************************************************/
/*                  calc_hann_wang_x0                                           */
/********************************************************************************/
/*  Description :   calculate the hand & wang initial value for mort trajectory */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   January 4, 2000                                             */
/*  Returns     :   void                                                        */
/*  Comments    :                                                               */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double          sdi             - stand density index                       */
/*  double          sdimx           - maximum stand density index               */
/*  double          rd_before       - relative density in last period           */
/*  double          rd_after        - relative density now                      */
/*  double          *x0             - the Hann & Wang x0 value                  */
/********************************************************************************/
/*  Formula :                                                                   */
/*  Source  : SYSTUM-1 + Hann & Wang Mortality paper                            */
/*  Coeffs  : Hann & Wang                                                       */
/********************************************************************************/
void calc_hann_wang_x0(   
   unsigned long   *return_code,
   double          sdi,
   double          bhtpa,
/*    double          qmd,  */
   double          sdimx,
   double			rd_before,
   double			rd_after,
   double          *x0 )
{

   double          ratio;

   /*  this should never happen*/
   if( bhtpa <= 0.0)
   {
      *x0 = 0.0;	
      return;
   }

   /* this should never happen */
   if( sdimx <= 0.0 )
   {
      *return_code = CONIFERS_MORT_ERROR; 
      *x0 = 0.0;
      return;
   }
	
   if (rd_after < CRITICAL_RD )
   {
      /* set to 0.0, no mortality and is currently < threshold */ 
      *x0 = 0.0;      
      *return_code = CONIFERS_SUCCESS; 
      return;         
   }


/*---------------------------------------------------------*/
/*   this is if it is now and was before above CRITICAL RD */
/*---------------------------------------------------------*/
   if (rd_before >= CRITICAL_RD && rd_after >= CRITICAL_RD)
   {
      /* x0 should be >0 here */
      if ( *x0 <= 0.0 )
      {
	 *return_code = CONIFERS_MORT_ERROR;
      }
      return; /* don't do anything leave x0 as it is */
   }



    
   if ( rd_after >= CRITICAL_RD && rd_before < CRITICAL_RD )
   {	
      /* this should never happen   */	
      if(*x0 > 0.0)
      {
	 *return_code = CONIFERS_MORT_ERROR;
	 return;
      }
      ratio = log(CRITICAL_RD) / log( sdi / sdimx );
      *x0   = log( bhtpa ) + 0.678688502f * log( ratio ); 
      /* for stands above the max sdi... */
      if( sdi >= sdimx ) 
      {
	 *x0 =    log( bhtpa ) +
	    0.678688502f * log(log(CRITICAL_RD) / log( 0.99f ) );
      }
   }
} 

/********************************************************************************/
/*                  calc_init_x0                                                */
/********************************************************************************/
/*  Description :   initialize x0 for initial read of stands>0.6mort trajectory */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   August 7, 2006                                              */
/*  Returns     :   void                                                        */
/*  Comments    :                                                               */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double          sdi             - stand density index                       */
/*  double          rd_current      - current rel density                       */
/*  double          bhtpa           - current bh trees per acre                 */
/*  double          sdimx           - maximum stand density index               */
/*  double          *x0             - the Hann & Wang x0 value                  */
/********************************************************************************/
/*  Formula :                                                                   */
/*  Source  : SYSTUM-1 + Hann & Wang Mortality paper                            */
/*  Coeffs  : Hann & Wang                                                       */
/********************************************************************************/
void calc_init_x0(   
   unsigned long   *return_code,
   double          sdi,
   double          rd_current,
   double          bhtpa,
   double          sdimx,
   double          *x0 )
{

   double          ratio;

   *return_code = CONIFERS_SUCCESS; 


   /*  this should never happen*/
   if( bhtpa <= 0.0 || sdi <= 0.0)
   {
      *x0 = 0.0;	
      return;
   }

   /* this should never happen */
   if( sdimx <= 0.0 )
   {
      *return_code = CONIFERS_MORT_ERROR; 
      *x0 = 0.0;
      return;
   }

   if (rd_current >= CRITICAL_RD)
   {
    
      ratio = log(CRITICAL_RD) / log( sdi / sdimx );
      *x0   = log( bhtpa ) + 0.678688502f * log( ratio ); 
      /* for stands above the max sdi... */
      if( sdi >= sdimx ) 
      {
	 *x0 =    log( bhtpa ) +
	    0.678688502f * log(log(CRITICAL_RD) / log( 0.99f ) );
      }
      return;
   }
} 

