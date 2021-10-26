#include "system_includes.h"
#include "state_struct.h"
#include "unique_ratio.h"
#include "get_ipsteps.h"
#include "output_combo_tsv.h"
#ifdef TM30
#include "get_steps.h"
#include "ies_tm_30_18.h"
#include "output_tm_30_results.h"
#include "output_tm_30_debug.h"
#endif
#ifdef DBG
#include "step_compare.h"
#endif
#include "process_combo.h"
void process_combo(struct state_struct *state,
		   struct loop_state_struct *loop_state, 
		   int nzc,
		   int64_t *num_in_tca_p,
		   int *istep_index,
		   double x,
		   double y,
		   double z,
		   double sqdist) {
  /*
    Determine if a set of power steps in istep_index is unique, 
    if (unique) {
      increment *num_in_tca_p
      if (state->select_only)
         if (state->dump_combos) {
	    print the power steps (istep_index(0:num_channels-1)
	    and the corresponding tristimulus values (x,y,z) 
            and the distance from the target center sqrt(sqdist)
	    to the state->dfp stream.
	 }
       } else {
         call ies_tm_30_18 to characterize the combination,
	 and update the summary data.
	 if (state->full_output) {
	    output the characterization data to the
	    state->ofp stream
         }
       }
     }

     Called by: eprocess_combos, bf_process_combos.
     Calls:     unique_ratio, get_ipsteps,output_combo_tsv,
                get_steps, ies_tm_30_18, output_tm_30_results,
		output_tm_30_debug, step_compare

     Uses the following fields of state:
        select_only,
	dump_combos,
	full_output,
	extra_log_info,
	dfp,
	ofp,
	lfp,
	rf_thresh,
	*source_rows,
	*ybar_t_channels,
	*e_t_channels,
	*sz_t_channels
	tm_30_results

     Sets the following field of state
        print_efit
	*tm_30_results

     Uses the following fields of loop_state:
        num_channels,
	max_index,
	*imap_sort,
	*unique_prime_factors,
	*upf_index,
	*fracts,
	steps,
	dbg_steps,
	ipstep_index,
	mm_istep_index

     Sets the following fields of loop_state
        *steps,
	*dbg_steps,
	*ipstep_index,
	*mm_istep_index,
	
  */

  double *tm_30_results;
  double *steps;
  double *fracts;
  double *ybar_t_channels;
  double *e_t_channels;
  double *sz_t_channels;

#ifdef DBG  
  double *dbg_steps;
  double *stest;
  double *channels;
#endif
  int *ipstep_index;
  int *mm_istep_index;
  int *imap_sort;
  int *unique_prime_factors;
  int *upf_index;

  double cct_vals[4];
  double rf_thresh;
  double ybar_dot_stest;
  double stest_sum;
  double ler;
  double sum_m;
  double sum_p;
  double m_over_p;
  double dist;

#ifdef DBG
  double thresh;
#endif

  int64_t num_in_tca;

  int dump_combos;
  int select_only;

  int num_channels;
  int max_index;

  int unique;
  int full_output;

  int i;
  int extra_log_info;

#ifdef DBG
  int    source_rows;
  int    ipad;

  int    j;
  int    i_dbg;
#endif  

  FILE *dfp;
  FILE *lfp;

  select_only          = state->select_only;
  dump_combos          = state->dump_combos;
  extra_log_info       = state->extra_log_info;
  dfp                  = state->dfp;
  lfp                  = state->lfp;


  num_channels 	       = loop_state->num_channels;
  max_index            = loop_state->max_index;

  ipstep_index         = loop_state->ipstep_index;
  mm_istep_index       = loop_state->mm_istep_index;
  imap_sort            = loop_state->imap_sort;
  unique_prime_factors = loop_state->unique_prime_factors;
  upf_index            = loop_state->upf_index;

  num_in_tca = *num_in_tca_p;

  /*
    Unique_ratio determines if we 
    have seen this ratio of istep_index values before.
    returning a unique flag = 1 if not,
    or 0 if we have. Also unique_ratio will,
    if the ratio is unique, convert it to its
    largest integer multiple in the power levels in
    mm_istep_index (for max-multiple)
  */
  unique = unique_ratio(num_channels,max_index,istep_index,
			unique_prime_factors,upf_index,mm_istep_index);

  if (unique) {
    num_in_tca += (int64_t)1;
    *num_in_tca_p = num_in_tca;
    if (select_only != 0) {
      if (dump_combos) {
	get_ipsteps(num_channels,imap_sort,mm_istep_index,ipstep_index);
	output_combo_tsv(dfp,num_channels,ipstep_index,x,y,z);
      
      }
    } else {
#ifdef TM30      
      full_output          = state->full_output;
      rf_thresh            = state->rf_thresh;
      ybar_t_channels      = state->ybar_t_channels;
      e_t_channels         = state->e_t_channels;
      sz_t_channels        = state->sz_t_channels;
      steps                = loop_state->step_fracts;
      fracts               = loop_state->fracts; /* set in init_loop_state */
      /* Now ies_tm_30_18 needs the steps in the unpermuted order */
      get_steps(num_channels,imap_sort,mm_istep_index,fracts,steps);

#ifdef DBG
      dbg_steps         = state->dbg_steps;
      state->print_efit = 0;
      source_rows       = state->source_rows;
      /*
	Right here we want to test steps for equality to dbg_steps
	and if so execute a print statemnt that allows us to stop on it.
	We only want to do this in debug mode.
	But we want the permuted steps.
      */
      thresh = 0.001;
      if (step_compare(num_channels,thresh,dbg_steps,steps) == 1) {
	fprintf(lfp,"Debug check");
	state->print_efit = 1;
      }
      if (state->print_efit) {
	stest = state->stest;
	for (i=0;i<source_rows;i++) {
	  stest[i] = 0.0;
	}
	channels = state->channels;
	for (j = 0;j<num_channels;j++) {
	  for (i=0;i<source_rows;i++) {
	    stest[i] += steps[j]*channels[i];
	  }
	  channels = channels + source_rows; /* Caution address arithmetic */
	}
	fprintf(lfp,"stest = \n");
	for (i=0;i<source_rows;i++) {
	  fprintf(lfp,"stest[%d] = %le\n",
		  i,stest[i]);
	}
	fflush(lfp);
      }
#endif
      tm_30_results       = state->tm_30_results;
      ies_tm_30_18(state,steps,tm_30_results,cct_vals);
      if (tm_30_results[0] >= rf_thresh) {
	if (full_output) {
	  ybar_dot_stest = ybar_t_channels[0] * steps[0];
	  for(i=1;i<num_channels;i++) {
	    ybar_dot_stest += ybar_t_channels[i] * steps[i];
	  }
	  stest_sum      = e_t_channels[0] * steps[0];
	  for (i=1;i<num_channels;i++) {
	    stest_sum    += e_t_channels[i] * steps[i];
	  }
	  ler   = 683.0*ybar_dot_stest/stest_sum;
	  sum_m = sz_t_channels[0] * steps[0];
	  for (i=1;i<num_channels;i++) {
	    sum_m += sz_t_channels[i] * steps[i];
	  }
	  sum_p = ybar_dot_stest;
	  m_over_p = sum_m/sum_p;
	  dist = sqrt(sqdist);
	  output_tm_30_results(state,loop_state,
			       steps,tm_30_results,ler,m_over_p,sum_p,dist,
			       cct_vals);
	  if (extra_log_info) {
	    get_ipsteps(num_channels,imap_sort,mm_istep_index,ipstep_index);
	    output_tm_30_debug(state,lfp,num_channels,ipstep_index,
			       dist,tm_30_results,nzc,ler,m_over_p,sum_p);
	  }
	} /* end if (full_output) */
      } /* end if (rf > rf_thresh) */
#endif /* TM30 code */     
    } /* end else not a select_only run */
  } /* end if (unique) */
}
