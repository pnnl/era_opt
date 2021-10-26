#include "system_includes.h"
#include "state_struct.h"
#include "groups_struct.h"
#include "comm_map_struct.h"
#include "ivtbar.h"
#include "ivtbxr.h"
#include "lb_summary.h"
int lb_summary(int nproc, struct state_struct *state) {
  /*
    Perform reductions on the num_in_tca, ix_delta_sum, iy_delta_sum fields
    of state, and the max_in_tca, max_iy_delta nad max_ix_delta fields
    These are the load_balance fields of state (hence the lb_).

    Called by: era_opt
    Calls:     ivtbar, ivtbxr

    Arguments:

    T = type   i = int, d = double, s = struct, f = file
    M = mode   0 = scalar, 1 = vector, 2 = array, * = pointer
    F = flow   i = input, o = output, b = input & output, w = scratch

    Variable   			Class  Description
               			TMF
   
    nproc                       i0i    number of processors.
    
    state                       s*b    state vector 
                    
    Fields of state used:
       num_in_tca;
       ix_delta_sum;
       iy_delta_sum;
       reduction_map
       lfp

    Fields of state set:
       ix_delta_sum,
       iy_delta_sum,
       num_in_tca,
       max_in_tca,
       max_iy_delta,
       max_ix_delta,
  */
  struct comm_map_struct *reduction_map;
  struct groups_struct *groups_state;
  int64_t *select_counts;
  int64_t *left_counts;
  int64_t *right_counts;
  int64_t *select_max;
  int64_t *left_max;
  int64_t *right_max;
  

  int number_counts;
  int i;

  int number_max;
  int select_only;
  
  int itag;
  int success;
  FILE *lfp;

  success          = 1;
  select_only      = state->select_only;
  groups_state     = state->groups_state;
  reduction_map    = state->reduction_map;
  lfp              = state->lfp;
  select_counts    = state->select_counts;
  left_counts      = state->ls_counts;
  right_counts     = state->rs_counts;
  select_max       = state->select_max;
  left_max         = state->ls_max;
  right_max        = state->rs_max;
  select_counts[0] = state->num_in_tca;
  select_counts[1] = state->iy_delta_sum;
  select_counts[2] = state->ix_delta_sum;
  select_counts[3] = 0;
  left_counts[0]   = (int64_t)0;
  left_counts[1]   = (int64_t)0;
  left_counts[2]   = (int64_t)0;
  left_counts[3]   = (int64_t)0;
  right_counts[0]  = (int64_t)0;
  right_counts[1]  = (int64_t)0;
  right_counts[2]  = (int64_t)0;
  right_counts[3]  = (int64_t)0;
  select_max[0]    = state->num_in_tca;
  /*
    NO I really want the maximum of their totals 
  select_max[1]    = state->max_iy_delta;
  select_max[2]    = state->max_ix_delta;
  */
  select_max[1]    = state->iy_delta_sum;
  select_max[2]    = state->ix_delta_sum;
  if (select_only == 0) {
    select_max[3]  = groups_state->number_results;
  } else {
    select_max[3]    = 0;
  }
  left_max[0]      = 0;
  left_max[1]      = 0;
  left_max[2]      = 0;
  left_max[3]      = 0;
  right_max[0]      = 0;
  right_max[1]      = 0;
  right_max[2]      = 0;
  right_max[3]      = 0;

  number_max       = 4;
  number_counts    = 4;
  /*
    Place neutral group counts and total number results into the all counts vector
    for doing the reduction sum accross processors.
  */
  if (nproc > 1) {
    itag = 701;
    success = ivtbar(number_counts,select_counts,left_counts,right_counts,itag,reduction_map,lfp);
    if (success != 1) {
      if (lfp) {
	fprintf(lfp,"lb_summary, zero return code for select_counts_reduction\n");
	fflush(lfp);
      }
    }
  } /* end if (nproc > 1) */
  if (success) {
    if (nproc > 1) {
      itag = 801;
      success = ivtbxr(number_max,select_max,left_max,right_max,itag,reduction_map,lfp);
      if (success != 1) {
	if (lfp) {
	  fprintf(lfp,"lb_summary, zero return code for select_max_reduction\n");
	  fflush(lfp);
	}
      }
    }
  }
  state->num_in_tca = select_counts[0];
  state->iy_delta_sum = select_counts[1];
  state->ix_delta_sum = select_counts[2];
  state->max_in_tca   = select_max[0];
  state->max_iy_delta_sum = select_max[1];
  state->max_ix_delta_sum = select_max[2];
  if (select_only == 0) {
    state->max_number_results = select_max[3];
  }
  return(success);
}
