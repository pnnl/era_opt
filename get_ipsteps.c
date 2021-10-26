#include "system_includes.h"
#include "state_struct.h"
#include "get_ipsteps.h"
void get_ipsteps(int num_channels, int *imap_sort, int *istep_index,
		 int *ipstep_index) {
  /*
struct loop_state_struct *loop_state) {
  */
  /*
    istep_index has a step set that lands in the target circle.
    multiply by the fracts and unpermute, for inptut to ies_tm_30_18
    and to print in the input order.
  */
  int i;
  int ipad;
  /*
  istep_index  = loop_state->istep_index;
  imap_sort    = loop_state->imap_sort;
  fracts       = loop_state->fracts;
  step_fracts  = loop_state->step_fracts;
  num_channels = loop_state->num_channels;
  */
  for (i=0;i<num_channels;i++) {
    ipstep_index[imap_sort[i]] = istep_index[i];
  }
}
