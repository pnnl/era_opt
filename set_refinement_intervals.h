#ifndef _SET_REFINEMENT_INTERVALS_H_
#define _SET_REFINEMENT_INTERVALS_H_ 1
extern void set_refinement_intervals(int number_scales,
				     int num_channels,
				     double *scales,
				     double *min_scale_p,
				     int *refinement_intervals, 
				     int *max_index_p,
				     int *istep_total_p, 
				     FILE *lfp);
#endif
