#include "system_includes.h"
#include "state_struct.h"
#include "set_refinement_intervals.h"
void set_refinement_intervals(int number_scales, int num_channels,
			      double *scales,
			      double *min_scale_p,
			      int *refinement_intervals, int *max_index_p,
			      int *istep_total_p, FILE *lfp) {
  /*
    Check to see that the scales are all multiples of the smallest scale,
    and set the reinfement_intervals array,the max_index value and istep_total,
    If the scales are not even multiples of the smallest scale set 
    istep_total to 0.
    It has already been verified in record scales that they are decreasing.
    Set istep total to be (max_step + 1) * num_channels where max_step
    is 1/smallest_scale.

    Called by: init_loop_state
  */
  double min_scale;
  double sratio;
  double extra;
  int max_index;
  int istep_total;
  int i;
  int success;
  success   = 1;
  min_scale = scales[number_scales-1];
  max_index = (int)(1.000001/min_scale);
  refinement_intervals[number_scales-1] = 1;
  if (number_scales > 1) {
    for (i=0;((i<=number_scales-2) && success);i++) {
      sratio = scales[i]/min_scale;
      refinement_intervals[i] = (int)(sratio + .00001);
      extra  = sratio - (double)refinement_intervals[i];
      if (extra < 0) extra = - extra;
      if (extra > .0000001) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"set_refinement_intervals: scale %le is not an integer"
		  " multiple of the smallest scale, %le\n", scales[i],min_scale);
	  fflush(lfp);
	}
      }
    }
  }
  if (success) {
    istep_total = (max_index + 1) * num_channels;
  } else {
    istep_total = 0;
  }
  *max_index_p   = max_index;
  *istep_total_p = istep_total;
  *min_scale_p   = min_scale;
  return;
}
