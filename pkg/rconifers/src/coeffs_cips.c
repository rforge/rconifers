/****************************************************************************/
/*                                                                          */
/*  cips_coeffs.c                                                           */
/*  coeffficients used in the CIPS variant of the CONIFERS growth model     */
/*                                                                          */
/****************************************************************************/

/* 	$Id: coeffs_cips.c 893 2012-06-13 17:41:46Z hamannj $ */

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"

/* todo: add/update version numbers for your model's coefficients   */
#define N_CIPS_CONIFERS_COEFFS      3
#define CIPS_COEFFS_VERSION      1.01
#define CIPS_MODEL_VERSION       1.01

/* commit test from new repos */

/* this structure MUST be the EXACT same structure def as the               */
/* one defined in the header file                                           */
struct CON_CIPS_COEFFS_RECORD 
{
    unsigned long   idx;                    /* functional species code      */
    unsigned long   type;                   /* 0, 1, 2, 3 or 4              */
    char            group[SP_LENGTH];       /* required by GUI, deprecated  */

    double  d6_growth[MAX_COEFFS];          /* 1  D1: growth coeffs         */
    double  ht_growth[MAX_COEFFS];          /* 2  D3: growth coeffs         */ 
    double  cr_growth[MAX_COEFFS];          /* 3  D4: hcb change            */
    double  crown_ratio[MAX_COEFFS];        /* 4  S3: crown ratio coeffs    */
    double  crown_width[MAX_COEFFS];        /* 5  S1: crown width coeffs    */
    double  max_crown_width[MAX_COEFFS];    /* 6  S2: max crown width       */
    double  d6_ht_dbh[MAX_COEFFS];          /* 7  S5, S7                    */
    double  d6_ht[MAX_COEFFS];              /* 8  S4, S9                    */
    double  dbh_ht[MAX_COEFFS];             /* 9  S6, S10                   */
    double  n_stems_from_ht[MAX_COEFFS];    /* 10 S11                       */
    double  dbh_growth[MAX_COEFFS];         /* 11 D5: dbh growth            */
    double  cfvolume4[MAX_COEFFS];          /* 12 S12 volume                */
    double  biomass[MAX_COEFFS];            /* 13 S13 biomass               */
    double  cw_growth[MAX_COEFFS];          /* 14 D6: cw_growth             */
    /* coefficients added for the CIPS variant */
    double	 dbh_ht_veg_cov[MAX_COEFFS];	/* dbh-ht prediction            */
    double	 dob_hi[MAX_COEFFS];            /* dob at height hi             */	
    double   mortality[MAX_COEFFS];         /* species specific mortality   */

} static cips_temp_coeffs[N_CIPS_CONIFERS_COEFFS] = {

/*                         0         1         2          3         4         5         6        7         8          9         10        11       12       13        14       15        */
{ /* index 0, FBR type = 2 , this is the default for brush */
0, SHRUB, "FBR",

/* these are the coeffs for the brush model */
/* d6_growth        */  { 0.052777, 0.046057,      0.0, 0.091103, 0.004174,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* ht_growth        */  { 0.376316,-3.461019, 1.677419,-0.099159,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* cr_growth        */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* crown_ratio      */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* crown_width      */  { 0.702099, 0.617866, 1.325870,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* max_crown_width  */  {      0.0,      0.0,     0.0,       0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* d6_ht_dbh        */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* d6_ht            */  {-0.958757, 0.536006,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* dbh_ht           */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* n_stems_from_ht  */  { 0.528185, 0.891697,-0.067543,-0.121772,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* dbh_growth       */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* cfvolume4        */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* biomass          */  {-6.909900, 2.854200, 0.078000,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* cw_growth        */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, 125.8000,  -1.1124, 0.074980, -8.32510, -0.06569},

/* dbh_ht_veg_cov   */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* dob_hi           */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* mortality        */  {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
},

/*                           0          1         2          3         4         5         6         7         8          9         10        11       12       13        14       15   */
{ /* index 1, FDF type = 0 */
1, CONIFER, "FDF",
/* d6_growth        */  {      0.0,       0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* ht_growth  2014  */  { 12.14808, -3.190740, 0.001047, -1.80396, -1.91691, -0.01265,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* cr_growth        */  { 0.170720,  0.080100, 0.259900,  0.31977, 0.085170, 3.465700, -1.78790, -0.69680, -0.64460, -0.64460,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* crown_ratio      */  { 6.069400, -0.000148, -1.62100, -0.10560,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* crown_width      */  { 0.335198,  0.494185, 1.227520,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* max_crown_width  */  { 4.636600,  1.607800,-0.009625,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* d6_ht_dbh        */  { 0.778122,  0.666671, 0.174612,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* d6_ht            */  {-1.375748,  1.016819,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* dbh_ht           */  { 7.127600, -5.364200,-0.261749,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* n_stems_from_ht  */  {     0.00,       0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* dbh_growth  2014 */  { 0.072151,  0.105527, -0.00258, -2.53834, -0.00451, 1.221857, -0.04529, 1.207514,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* cfvolume4        */  { 0.002244,  1.943420, 0.996400, 0.999200,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* biomass          */  {      0.0,       0.0,      0.0,      0.0,      0.0,      0.0, -2.84620,  1.70090, -3.69410,  2.13820, -3.03960, 2.595100,      0.0,      0.0,     0.0,      0.0},
/* cw_growth (1)    */  { 0.860817,  0.957470, -0.16419,      0.0,-0.098730,      0.0,      0.0, 0.100000,      0.0,      0.0,      0.0, 125.8000,  -1.1124,  0.07498, -8.3251, -0.06569},

/* dbh_ht_veg_cov   */  { -10.4003,   2.44412, -0.04612, 0.061653,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* dob_hi           */  {   0.4361,   -0.0504,   0.8046, -11.8347,   5.7161,   0.2206,  11.3148,   1.9433,   6.3261,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},
/* mortality        */  {  -2.8129,   -0.4986,   0.9901,  -0.6312,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,      0.0},

},

/*                         0         1         2          3         4         5         6        7         8          9         10        11       12       13        14       15        */
{ /* index 16, FNS type = 4 */
2, NON_STOCKED, "FNS",
/* d6_growth        */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* ht_growth        */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* cr_growth        */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* crown_ratio      */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* crown_width      */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* max_crown_width  */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* d6_ht_dbh        */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* d6_ht            */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* dbh_ht           */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* n_stems_from_ht  */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* dbh_growth       */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* cfvolume4        */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* biomass          */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* cw_growth (1)    */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},

/* dbh_ht_veg_cov   */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* dob_hi           */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
/* mortality        */{      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0},
}

};



/* This function allocates an array of n_coeffs records     */
/* which is the array that will be returned. The function   */
/* then copies the local coeffs into that array and the     */
/* returns                                                  */
struct COEFFS_RECORD *con_cips_init_coeffs( 
   unsigned long               *n_coeffs, 
   double                       *coeffs_version,
   double                       *model_version )
{

   struct COEFFS_RECORD *temp_coeffs; 
   
   /* todo: make sure the size of the coefficient records match */
    if( sizeof( struct COEFFS_RECORD ) == sizeof( struct CON_CIPS_COEFFS_RECORD ) )
    {
        /* these are defined by the model and coefficients version */
        *n_coeffs       = N_CIPS_CONIFERS_COEFFS;
        *coeffs_version = (double)CIPS_COEFFS_VERSION;
        *model_version  = (double)CIPS_MODEL_VERSION;

        /* allocate the memory to hold the coefficients */
        temp_coeffs = (struct COEFFS_RECORD *)calloc( 
        *n_coeffs, sizeof( struct COEFFS_RECORD ) );

        /* and copy the coefficients into the return array */
        memcpy( temp_coeffs, cips_temp_coeffs, sizeof( cips_temp_coeffs ) );
    }
    else
    {
        *n_coeffs       = 0;
        *coeffs_version = 0.0;
        *model_version  = 0.0;
        return NULL;
    }

   return temp_coeffs;
}

