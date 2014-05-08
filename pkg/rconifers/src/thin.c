/********************************************************************************/
/*                                                                              */
/*  thin.c                                                                      */
/*  functions used to do thin and apply mortality to the tree/plant list        */
/*  since mortality (sdi driven) is applied in much the same way as a thinning  */
/*  this is the logical place to keep this stuff                                */
/*                                                                              */
/********************************************************************************/

/* 	$Id: thin.c 930 2014-01-29 21:48:28Z mritchie $	 */

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"


static int compare_DBH_increasing( 
    const void *ptr1, 
    const void *ptr2 );

/* these functions have been made static since they don't need to be exposed    */
/* the thin_plot function is the external interface                             */
static void do_expf_sp_thinning( 
    unsigned long           *return_code,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,
    struct COEFFS_RECORD    *coeffs_ptr,
    unsigned long           thin_species_idx,
    double                  target,
    double                  *plants_removed, 
    double                  *ba_removed);

static void do_expf_thinning( 
    unsigned long           *return_code,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  target,
    double                  *plants_removed, 
    double                  *ba_removed);

static void do_expf_thinning_from_below( 
    unsigned long           *return_code,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  target,
    double                  *plants_removed, 
    double                  *ba_removed);

static void do_expf_sp_thinning_from_below( 
    unsigned long           *return_code,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,
    struct COEFFS_RECORD    *coeffs_ptr,
    unsigned long           thin_species_idx,
    double                  target,
    double                  *plants_removed, 
    double                  *ba_removed);

/****************************************************************************/
/* sorting function for tree records                                        */
/****************************************************************************/
static int compare_DBH_increasing( 
    const void *ptr1, 
    const void *ptr2 )
{
    struct PLANT_RECORD *p1_ptr;
    struct PLANT_RECORD *p2_ptr;

    p1_ptr = (struct PLANT_RECORD*)ptr1;
    p2_ptr = (struct PLANT_RECORD*)ptr2;

	if( p1_ptr->dbh < p2_ptr->dbh )
    {
        return -1;
    }
	if( p1_ptr->dbh > p2_ptr->dbh )
    {
        return 1;
    }
	else
    {
        return 0;
    }

}


/********************************************************************************/
/* thin_plot                                                                    */
/********************************************************************************/
/*  Description :   thins the plants on the plot either for mortality or        */
/*                  management                                                  */
/*  Author      :   Martin Ritchie                                              */
/*  Date        :   January 14, 2000                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   This is the main thinning function                          */
/*                   it is used to thin plants out according to provided        */
/*                   instructions or to apply sdi related mortality             */
/*                   it is not used for individual-plant or endemic mortality,  */
/*                   as that happens as a part of the growth process            */
/*                   it is called in a loop which runs throught the plots, so   */
/*                   when you get here, you already know what plot you are on   */
/*  Arguments   :                                                               */
/*     unsigned long         *return_code   - pointer to a return code          */
/*     unsigned long         n_plants       - total number fo plants in the     */
/*                                            plants pointer array              */
/*     struct PLANT_RECORD   *plants_ptr    - array of plants in the  sample to */
/*                                            sample to be projected            */
/*     struct PLOT_RECORD    *plot_ptr      - pointer to the current plot that  */
/*                                               is to be grown                 */
/*     unsigned long         n_species      - size of the species_ptr           */
/*     struct SPECIES_RECORD *species_ptr   - array of SPECIES_RECORD's that    */
/*                                            hold species specific information */
/*     unsigned long         target         - species indicator                 */
/*     int                   thin_guide     - tells how to take plants:         */
/*                                              proportional, below, above, etc */
/*     double                target         - tells how much to remove, or leave*/
/*                                                                              */
/*                                                                              */
/********************************************************************************/
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
    double                  *ba_removed )
{

    /* the variable: thin_species_code is the species to be thinned,  */
    /* it is NULL for all species or sdi driven mortality             */

    unsigned long   i;
    struct  PLANT_RECORD    *plant_ptr;
    struct COEFFS_RECORD    *c_ptr;

    /*  declare variables  */
    //unsigned long           start_idx; removed, unused,  by mwr jan 2014
    //unsigned long           end_idx;   removed, unused,  by mwr jan 2014
    //unsigned long           n_plant_records_on_plot; rem, unused,  mwr jan 2014
    //double                  thinning_proportion;rem, unused,  mwr jan 2014
    //double                  plants_removed;
    //double                  ba_removed;
    //double                  total_trees;

    /*  MOD001 initialize variables    */
    // next two variables unused, commented out Jan 2014
    //start_idx                = 0;   /*  the starting index for current plot */
    //end_idx                  = 0;   /*  the ending index for current plot   */
    //n_plant_records_on_plot  = 0;   /*  number of plant records on plot     */
    //thinning_proportion      = 0.0; /*  the actual ratio to be removed      */
    //total_trees              = 0.0;  /* total trees per acre on the plot   */

    *plants_removed           = 0.0;
    *ba_removed               = 0.0;

    switch( thin_guide )
    {
        /* this is the routine that actually does the sdi mortality */
        case DO_SDI_MORT:    

            plant_ptr = &plants_ptr[0];
            for( i = 0; i < n_plants; i++, plant_ptr++ )
            {
                /* if you're not on the plot, don't thin it */
                if( plant_ptr->plot != plot_ptr->plot )
                {
                    continue;
                }
                c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

                if( is_tree( c_ptr ) )
                {       
                    /* target is mortality proportion from sdi mortality    */
                    plant_ptr->expf_change  = plant_ptr->expf * target;
                
                    plant_ptr->expf         -= plant_ptr->expf_change;
                }
            }
        break;

        case DO_EXPF_SP_THIN:                    
            /* TPA thin for a species             */
            do_expf_sp_thinning( 
                            return_code,
                            n_plants,
                            plants_ptr,
                            plot_ptr,
                            n_species,
                            species_ptr,
                            n_coeffs,
                            coeffs_ptr,
                            //thin_species_code,
                            thin_species_idx,
                            target,
                            plants_removed, /* effectively returns number thinned */
                            ba_removed);

        break;

        case DO_EXPF_THIN:        
            /*  thin to target tpa for all species on the plot  */
            do_expf_thinning(
                            return_code,
                            n_plants,
                            plants_ptr,
                            plot_ptr,
                            n_species,
                            species_ptr,
                            n_coeffs,
                            coeffs_ptr,
                            target,
                            plants_removed, 
                            ba_removed );
        break;     


        /*  thin to target tpa for all species on the plot  */
        case DO_EXPF_THIN_FROM_BELOW:
            do_expf_thinning_from_below(
                            return_code,
                            n_plants,
                            plants_ptr,
                            plot_ptr,
                            n_species,
                            species_ptr,
                            n_coeffs,
                            coeffs_ptr,
                            target,
                            plants_removed, 
                            ba_removed );

        break;     

        /*  thin to target tpa for specific species on the plot  */
        case DO_EXPF_SP_THIN_FROM_BELOW:
            do_expf_sp_thinning_from_below(
                            return_code,
                            n_plants,
                            plants_ptr,
                            plot_ptr,
                            n_species,
                            species_ptr,
                            n_coeffs,
                            coeffs_ptr,
                            //thin_species_code,
                            thin_species_idx,
                            target,
                            plants_removed, 
                            ba_removed );

        break;     


        /* todo: if you want to create a new thinning type, add it here. */


        default:
            *return_code=THINNING_ERROR;
        break;
    }
}

/*  MOD004  */
/********************************************************************************/
/* do_expf_sp_thinning                                                          */
/********************************************************************************/
/*  Description :   thins the plants on to a specified tpa for a specified      */
/*                   species                                                    */
/*  Author      :   Martin Ritchie                                              */
/*  Date        :   January 20, 2000                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   This is the thinning function for species tpa thinning      */
/*                   it is called only from thin_plot                           */
/*  Arguments   :                                                               */
/*     unsigned long         *return_code   - pointer to a return code          */
/*     unsigned long         n_plants       - total number fo plants in the     */
/*                                            plants pointer array              */
/*     struct PLANT_RECORD   *plants_ptr    - array of plants in the  sample to */
/*                                            sample to be projected            */
/*     struct PLOT_RECORD    *plot_ptr      - pointer to the current plot that  */
/*                                               is to be grown                 */
/*     unsigned long         n_species      - size of the species_ptr           */
/*     struct SPECIES_RECORD *species_ptr   - array of SPECIES_RECORD's that    */
/*                                            hold species specific information */
/*     unsigned long         n_coeffs         number of coefficients            */
/*     struct COEFFS_RECORD  *coeffs_ptr    - pointer to the coefficients       */
/*     char                  thin_species_code - species to be thinned          */
/*     double                target         - tells how much to remove, or leave*/
/*     double                *plants_removed - returns number of stems thinned  */
/*     double                *ba_removed     - returns basal area (d6) removed  */
/********************************************************************************/
static void do_expf_sp_thinning( 
    unsigned long           *return_code,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,
    struct COEFFS_RECORD    *coeffs_ptr,
    unsigned long           thin_species_idx,
    double                  target,
    double                  *plants_removed, 
    double                  *ba_removed)     
{

/* declarations  */
    unsigned long           start_idx;
    unsigned long           end_idx;
    unsigned long           n_plant_records_on_plot;
    double                  thinning_proportion;
    unsigned long           i;
    struct  PLANT_RECORD    *plant_ptr;

    unsigned long           n_sp_on_plot;       /* size of the species sums array */
    struct SUMMARY_RECORD   *plot_sum_ptr;      /* array of species sums */
    struct SUMMARY_RECORD   *target_species_plot_summary;


    /* initializations */
    start_idx               =0;   /* starting place in tree list for current plot */
    end_idx                 =0;   /* ending place in tree list for current plot   */
    n_plant_records_on_plot =0;   /* number of plants on the current plot         */
    thinning_proportion     =0.0; /* proportion left after thinning               */
    n_sp_on_plot            =0;   /* number of species on the current plot        */
    *plants_removed          =0.0;
    *ba_removed              =0.0;

    /*  This returns num spec on the plot */
    get_plant_indecies_for_plot( 
                            return_code,
                            plot_ptr,
                            n_plants,
                            plants_ptr,
                            &start_idx,      
                            &end_idx,
                            &n_plant_records_on_plot );

    /* build species summary array */
    plot_sum_ptr = build_species_summaries( 
                            return_code,
                            n_species,
                            species_ptr,
                            n_plant_records_on_plot,
                            &plants_ptr[start_idx],    
                            &n_sp_on_plot);            

    update_species_summaries( 
                            return_code,
                            n_species,
                            species_ptr,
                            n_coeffs,
                            coeffs_ptr,
                            n_plant_records_on_plot,
                            &plants_ptr[start_idx],
                            1,
                            n_sp_on_plot,
                            plot_sum_ptr );


    target_species_plot_summary = get_summary_from_code(    n_sp_on_plot,
                                                            plot_sum_ptr,
                                                            thin_species_idx );

    /*species not on plot */
    if(target_species_plot_summary == NULL)  
    /* target species not found here   */
    {

        *return_code=INVALID_SP_CODE;
        free(plot_sum_ptr);
        return;
    }

    /* then there are not enough plants here to thin  */
    if(target_species_plot_summary->expf <= target)
    {
        free(plot_sum_ptr);
        return;
    }

    thinning_proportion = target / target_species_plot_summary->expf;
    free(plot_sum_ptr);

    /*  now loop through and thin the trees on this plot  */
    plant_ptr = &plants_ptr[start_idx];
    for( i = start_idx; i <= end_idx; i++, plant_ptr++ )
    {
        /* if you're not on the plot, don't thin it */
        if( plant_ptr->plot != plot_ptr->plot )
        {
            continue;
        }

        if( plant_ptr->sp_idx != thin_species_idx )
        {
            continue;
        }

        /* target is mortality proportion from sdi mortality*/
        plant_ptr->expf_change  = 
                   plant_ptr->expf * (1.0 - thinning_proportion);
        plant_ptr->expf         *= thinning_proportion;
        /* tally the number and ba removed */
        *plants_removed          += plant_ptr->expf_change;
        *ba_removed              +=
                   plant_ptr->expf_change * plant_ptr->d6_area;
    }
}


/*  MOD005  */
/********************************************************************************/
/* do_expf_thinning                                                             */
/********************************************************************************/
/*  Description :   thins the plants on to a specified tpa for all species      */
/*                   of trees                                                   */
/*  Author      :   Martin Ritchie                                              */
/*  Date        :   January 20, 2000                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   This is the thinning function for species tpa thinning      */
/*                   it is called only from thin_plot                           */
/*  Arguments   :                                                               */
/*     unsigned long         *return_code   - pointer to a return code          */
/*     unsigned long         n_plants       - total number fo plants in the     */
/*                                            plants pointer array              */
/*     struct PLANT_RECORD   *plants_ptr    - array of plants in the  sample to */
/*                                            sample to be projected            */
/*     struct PLOT_RECORD    *plot_ptr      - pointer to the current plot that  */
/*                                               is to be grown                 */
/*     unsigned long         n_species      - size of the species_ptr           */
/*     struct SPECIES_RECORD *species_ptr   - array of SPECIES_RECORD's that    */
/*                                            hold species specific information */
/*     unsigned long         n_coeffs         number of coefficients            */
/*     struct COEFFS_RECORD  *coeffs_ptr    - pointer to the coefficients       */
/*     double                target         - tells how much to remove, or leave*/
/*     double                *plants_removed - returns number of stems thinned  */
/*     double                *ba_removed    - returns basal area (d6) removed   */
/********************************************************************************/
static void do_expf_thinning( 
    unsigned long           *return_code,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  target,
    double                  *plants_removed, 
    double                  *ba_removed)
{

    /* declarations  */
    unsigned long           start_idx;
    unsigned long           end_idx;
    unsigned long           n_plant_records_on_plot;
    double                  thinning_proportion;
    unsigned long           i;
    struct  PLANT_RECORD    *plant_ptr;


    struct COEFFS_RECORD    *c_ptr;           /* temporary coeffs record */
    unsigned long           n_sp_on_plot;     /* size of spec sums array */
    struct SUMMARY_RECORD   *plot_sum_ptr;    /* array of species sums   */
    struct  SUMMARY_RECORD  *temp_plot_record;
    double                  total_trees;

/* initializations */
    start_idx               =0;   /* starting place in tree list for plot */
    end_idx                 =0;   /* ending place in tree list for plot   */
    n_plant_records_on_plot =0;   /* number of plants on the plot         */
    thinning_proportion     =0.0; /* proportion left after thinning       */
    n_sp_on_plot            =0;   /* number of species on the plot        */
    *plants_removed         =0.0;
    *ba_removed             =0.0;
    total_trees             =0.0;

    /*  This returns num spec on the plot */
    get_plant_indecies_for_plot(
                            return_code,
                            plot_ptr,
                            n_plants,
                            plants_ptr,
                            &start_idx,
                            &end_idx,
                            &n_plant_records_on_plot );
    /* build species summary array */
    plot_sum_ptr = build_species_summaries( 
                            return_code,
                            n_species,
                            species_ptr,
                            n_plant_records_on_plot,
                            &plants_ptr[start_idx],    
                            &n_sp_on_plot);            

    update_species_summaries( 
                            return_code,
                            n_species,
                            species_ptr,
                            n_coeffs,
                            coeffs_ptr,
                            n_plant_records_on_plot,
                            &plants_ptr[start_idx],
                            1,
                            n_sp_on_plot,
                            plot_sum_ptr );

    temp_plot_record = &plot_sum_ptr[0];

    /*loop through the specs */
    for ( i=0; i < n_sp_on_plot; i++, temp_plot_record++ )
    {
        c_ptr = &coeffs_ptr[species_ptr[temp_plot_record->code].fsp_idx];

        /* check for trees only!  */
        if( is_tree( c_ptr ) )
        {
            total_trees += temp_plot_record->expf;
        }
    }

    /* then there are not enough plants here to thin  */
    if( total_trees <= target )
    {   
        free(plot_sum_ptr);
        return;
    }
    /* set thinning proportion */
    thinning_proportion = target / total_trees;
    free(plot_sum_ptr);
    plot_sum_ptr=NULL;

    plant_ptr = &plants_ptr[start_idx];
    for( i = start_idx; i <= end_idx; i++, plant_ptr++ )
    {
        c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

        /* check for trees only!  */
        if( is_tree( c_ptr ) )
        {
            plant_ptr->expf_change  = 
                    plant_ptr->expf * ( 1.0 - thinning_proportion );
            plant_ptr->expf         *= thinning_proportion;

            /* tally the number and ba removed */
            *plants_removed         += plant_ptr->expf_change;
            *ba_removed             +=
                    plant_ptr->expf_change * plant_ptr->d6_area;

            /* this next line is here for debugging purposes */
            total_trees             -= plant_ptr->expf_change;
        }
    }
}






/********************************************************************************/
/* do_expf_thinning_from_below                                                  */
/********************************************************************************/
/*  Description :   thins the plants on to a specified tpa for all species      */
/*                   of trees                                                   */
/*  Author      :   Martin Ritchie                                              */
/*  Date        :   January 20, 2000                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   This is the thinning function for species tpa thinning      */
/*                   it is called only from thin_plot                           */
/*  Arguments   :                                                               */
/*     unsigned long         *return_code   - pointer to a return code          */
/*     unsigned long         n_plants       - total number fo plants in the     */
/*                                            plants pointer array              */
/*     struct PLANT_RECORD   *plants_ptr    - array of plants in the  sample to */
/*                                            sample to be projected            */
/*     struct PLOT_RECORD    *plot_ptr      - pointer to the current plot that  */
/*                                               is to be grown                 */
/*     unsigned long         n_species      - size of the species_ptr           */
/*     struct SPECIES_RECORD *species_ptr   - array of SPECIES_RECORD's that    */
/*                                            hold species specific information */
/*     unsigned long         n_coeffs         number of coefficients            */
/*     struct COEFFS_RECORD  *coeffs_ptr    - pointer to the coefficients       */
/*     double                target         - tells how much to remove, or leave*/
/*     double                *plants_removed - returns number of stems thinned  */
/*     double                *ba_removed    - returns basal area (d6) removed   */
/********************************************************************************/
static void do_expf_thinning_from_below( 
    unsigned long           *return_code,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,
    struct COEFFS_RECORD    *coeffs_ptr,
    double                  target,
    double                  *plants_removed, 
    double                  *ba_removed)
{

    /* declarations  */
    unsigned long           start_idx;
    unsigned long           end_idx;
    unsigned long           n_plant_records_on_plot;
    double                  thinning_proportion;
    unsigned long           i;
    struct  PLANT_RECORD    *plant_ptr;


    struct COEFFS_RECORD    *c_ptr;           /* temporary coeffs record */
    unsigned long           n_sp_on_plot;     /* size of spec sums array */
    struct SUMMARY_RECORD   *plot_sum_ptr;    /* array of species sums   */
    struct  SUMMARY_RECORD  *temp_plot_record;
    double                  total_trees;
    double                  tree_value;
    double                  thin_value;

    /* initializations */
    start_idx               =0;   /* starting place in tree list for plot */
    end_idx                 =0;   /* ending place in tree list for plot   */
    n_plant_records_on_plot =0;   /* number of plants on the plot         */
    thinning_proportion     =0.0; /* proportion left after thinning       */
    n_sp_on_plot            =0;   /* number of species on the plot        */
    *plants_removed         =0.0;
    *ba_removed             =0.0;
    total_trees             =0.0;

    /*  This returns num spec on the plot */
    get_plant_indecies_for_plot(    return_code,
                                    plot_ptr,
                                    n_plants,
                                    plants_ptr,
                                    &start_idx,
                                    &end_idx,
                                    &n_plant_records_on_plot );

    /* build species summary array */
    plot_sum_ptr = build_species_summaries( return_code,
                                            n_species,
                                            species_ptr,
                                            n_plant_records_on_plot,
                                            &plants_ptr[start_idx],    
                                            &n_sp_on_plot);            

    /* compute the species summaries for the array */
    update_species_summaries(   return_code,
                                n_species,
                                species_ptr,
                                n_coeffs,
                                coeffs_ptr,
                                n_plant_records_on_plot,
                                &plants_ptr[start_idx],
                                1,
                                n_sp_on_plot,
                                plot_sum_ptr );


    /* loop through the specs, because you need to      */
    /* know the number of trees per acre on the plot    */
    temp_plot_record = &plot_sum_ptr[0];
    for( i = 0; i < n_sp_on_plot; i++, temp_plot_record++ )
    {
        c_ptr = &coeffs_ptr[species_ptr[temp_plot_record->code].fsp_idx];

        /* check for trees only!  */
        if( is_tree( c_ptr ) )
        {
            total_trees += temp_plot_record->expf;
        }
    }

    /* get rid of the plot summary poitner */
    free( plot_sum_ptr );
    plot_sum_ptr = NULL;

    /* then there are not enough plants here to thin  */
    if( total_trees <= target )
    {   
        return;
    }

    /* sort the plot by increasing diameters */
    qsort(  (void*)&plants_ptr[start_idx], 
            (size_t)(n_plant_records_on_plot), 
            sizeof( struct PLANT_RECORD ),
		    compare_DBH_increasing );    


    thin_value      = total_trees - target;

    /* set thinning proportion */
    //thinning_proportion = target / total_trees;
    //thin_tally = 0.0;
    plant_ptr = &plants_ptr[start_idx];
    for( i = start_idx; i <= end_idx; i++, plant_ptr++ )
    {
        c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

        /* check for trees only!  */
        if( is_tree( c_ptr ) )
        {
            /* new code */
            tree_value = plant_ptr->expf;
            
            /* modify the tree record record EXPF accordingly   */
            if( tree_value < ( thin_value - *plants_removed ) )
            {
                thinning_proportion     = 0.0;
                plant_ptr->expf_change  = 
                    plant_ptr->expf * ( 1.0 - thinning_proportion );
                plant_ptr->expf         *= thinning_proportion;

                *plants_removed         += plant_ptr->expf_change;
                *ba_removed             +=
                    plant_ptr->expf_change * plant_ptr->d6_area;
            }
            else
            {
                thinning_proportion     = ( thin_value - *plants_removed ) / tree_value;

                plant_ptr->expf_change  = 
                    plant_ptr->expf * ( 1.0 - thinning_proportion );
                plant_ptr->expf         *= thinning_proportion;

                *plants_removed         += plant_ptr->expf_change;
                *ba_removed             +=
                    plant_ptr->expf_change * plant_ptr->d6_area;

                break;
            }
        }
    }
}




/********************************************************************************/
/* do_expf_sp_thinning_from_below                                               */
/********************************************************************************/
/*  Description :   thins the plants on to a specified tpa for all species      */
/*                   of trees                                                   */
/*  Author      :   Martin Ritchie                                              */
/*  Date        :   January 20, 2000                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   This is the thinning function for species tpa thinning      */
/*                   it is called only from thin_plot                           */
/*  Arguments   :                                                               */
/*     unsigned long         *return_code   - pointer to a return code          */
/*     unsigned long         n_plants       - total number fo plants in the     */
/*                                            plants pointer array              */
/*     struct PLANT_RECORD   *plants_ptr    - array of plants in the  sample to */
/*                                            sample to be projected            */
/*     struct PLOT_RECORD    *plot_ptr      - pointer to the current plot that  */
/*                                               is to be grown                 */
/*     unsigned long         n_species      - size of the species_ptr           */
/*     struct SPECIES_RECORD *species_ptr   - array of SPECIES_RECORD's that    */
/*                                            hold species specific information */
/*     unsigned long         n_coeffs         number of coefficients            */
/*     struct COEFFS_RECORD  *coeffs_ptr    - pointer to the coefficients       */
/*     double                target         - tells how much to remove, or leave*/
/*     double                *plants_removed - returns number of stems thinned  */
/*     double                *ba_removed    - returns basal area (d6) removed   */
/********************************************************************************/
static void do_expf_sp_thinning_from_below( 
    unsigned long           *return_code,
    unsigned long           n_plants,
    struct PLANT_RECORD     *plants_ptr,
    struct PLOT_RECORD      *plot_ptr,
    unsigned long           n_species,
    struct SPECIES_RECORD   *species_ptr,
    unsigned long           n_coeffs,
    struct COEFFS_RECORD    *coeffs_ptr,
    unsigned long           thin_species_idx,
    double                  target,
    double                  *plants_removed, 
    double                  *ba_removed)
{

    /* declarations  */
    unsigned long           start_idx;
    unsigned long           end_idx;
    unsigned long           n_plant_records_on_plot;
    double                  thinning_proportion;
    unsigned long           i;
    struct  PLANT_RECORD    *plant_ptr;

    struct COEFFS_RECORD    *c_ptr;           /* temporary coeffs record */
    unsigned long           n_sp_on_plot;     /* size of spec sums array */
    struct SUMMARY_RECORD   *plot_sum_ptr;    /* array of species sums   */
    struct  SUMMARY_RECORD  *temp_plot_record;
    double                  total_trees;
    double                  tree_value;
    double                  thin_value;

    /* initializations */
    start_idx               =0;   /* starting place in tree list for plot */
    end_idx                 =0;   /* ending place in tree list for plot   */
    n_plant_records_on_plot =0;   /* number of plants on the plot         */
    thinning_proportion     =0.0; /* proportion left after thinning       */
    n_sp_on_plot            =0;   /* number of species on the plot        */
    *plants_removed         =0.0;
    *ba_removed             =0.0;
    total_trees             =0.0;

    /*  This returns num spec on the plot */
    get_plant_indecies_for_plot(    return_code,
                                    plot_ptr,
                                    n_plants,
                                    plants_ptr,
                                    &start_idx,
                                    &end_idx,
                                    &n_plant_records_on_plot );

    /* build species summary array */
    plot_sum_ptr = build_species_summaries( return_code,
                                            n_species,
                                            species_ptr,
                                            n_plant_records_on_plot,
                                            &plants_ptr[start_idx],    
                                            &n_sp_on_plot);            

    /* compute the species summaries for the array */
    update_species_summaries(   return_code,
                                n_species,
                                species_ptr,
                                n_coeffs,
                                coeffs_ptr,
                                n_plant_records_on_plot,
                                &plants_ptr[start_idx],
                                1,
                                n_sp_on_plot,
                                plot_sum_ptr );


    /* loop through the specs, because you need to      */
    /* know the number of trees per acre on the plot    */
    temp_plot_record = &plot_sum_ptr[0];
    for( i = 0; i < n_sp_on_plot; i++, temp_plot_record++ )
    {
        c_ptr = &coeffs_ptr[species_ptr[temp_plot_record->code].fsp_idx];

        if( is_tree( c_ptr ) && ( temp_plot_record->code == thin_species_idx ) )
        {
            total_trees += temp_plot_record->expf;
        }
    }

    /* get rid of the plot summary poitner */
    free( plot_sum_ptr );
    plot_sum_ptr = NULL;

    /* then there are not enough plants here to thin  */
    if( total_trees <= target )
    {   
        return;
    }

    /* sort the plot by increasing diameters */
    qsort(  (void*)&plants_ptr[start_idx], 
            (size_t)(n_plant_records_on_plot), 
            sizeof( struct PLANT_RECORD ),
		    compare_DBH_increasing );    


    thin_value      = total_trees - target;

    /* set thinning proportion */
    //thinning_proportion = target / total_trees;
    //thin_tally = 0.0;
    plant_ptr = &plants_ptr[start_idx];
    for( i = start_idx; i <= end_idx; i++, plant_ptr++ )
    {
        c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

        if( is_tree( c_ptr ) && ( plant_ptr->sp_idx == thin_species_idx ) )
        {

            /* new code */
            tree_value = plant_ptr->expf;
            
            /* modify the tree record record EXPF accordingly   */
            if( tree_value < ( thin_value - *plants_removed ) )
            {
                thinning_proportion     = 0.0;
                plant_ptr->expf_change  = 
                    plant_ptr->expf * ( 1.0 - thinning_proportion );
                plant_ptr->expf         *= thinning_proportion;

                *plants_removed         += plant_ptr->expf_change;
                *ba_removed             +=
                    plant_ptr->expf_change * plant_ptr->d6_area;
            }
            else
            {
                thinning_proportion     = ( thin_value - *plants_removed ) / tree_value;

                plant_ptr->expf_change  = 
                    plant_ptr->expf * thinning_proportion;
                plant_ptr->expf         *= ( 1.0 - thinning_proportion );

                *plants_removed         += plant_ptr->expf_change;
                *ba_removed             +=
                    plant_ptr->expf_change * plant_ptr->d6_area;

                break;
            }

        }

    }
}




