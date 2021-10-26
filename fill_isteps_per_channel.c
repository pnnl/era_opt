#include "system_includes.h"
#include "state_struct.h"
#include "fill_isteps_per_channel.h"
void fill_isteps_per_channel(struct state_struct *state) {
  /*
    Fill out an under specified isteps_per_channel array propagating
    the last-set value to the unsedt values.
    Called by: rthous_opt

    Uses the following fields of state
       num_channels
       nc_steps_set,
       isteps_per_channel[0:nc_steps_set-1]

    Sets the following fields of state
       isteps_per_channel[nc_steps_set:num_channels-1]
  */
  int *isteps_per_channel;
  int num_channels;
  int nc_steps_set;
  int last_steps_per_channel;
  int i;
  isteps_per_channel = state->isteps_per_channel;
  num_channels       = state->num_channels;
  nc_steps_set       = state->nc_steps_set;
  if (nc_steps_set < num_channels) {
    if (nc_steps_set > 0) {
      last_steps_per_channel = isteps_per_channel[nc_steps_set-1];
      /*
	if (lfp) {
	fprintf(lfp,"last_steps_per_channel = %d\n",last_steps_per_channel);
	fflush(lfp);
	}
      */
      for (i=nc_steps_set;i<num_channels;i++) {
	isteps_per_channel[i] = last_steps_per_channel;
      }
    } else {
      last_steps_per_channel = 6;
      for (i=0;i<num_channels;i++) {
	isteps_per_channel[i] = last_steps_per_channel;
      }
    }
  }
  return;
}
