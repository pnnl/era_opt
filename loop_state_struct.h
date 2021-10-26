#ifndef _LOOP_STATE_STRUCT_H_
#define _LOOP_STATE_STRUCT_H_ 1
struct loop_state_struct {
  void   *task_msg;
  double *stepped_channel_xyz; /* 3 * sum(steps_per_channel) */
  double *prev_xyz; /* 3 * num_channels */
  double *xyz_rps; /* 3 * num_channels , reverse partials sums of 
		      tristimulus xyz */
  double *curr_xyz; /* 3 */
  double *xyz_sort; /* 6 * num_channels */
  double *fracts;    /* num_channels  fracts[i] = 1/(msteps_per_channel[i]-1)
			for i = 0,num_channels-1 */
  double *pstep_fracts; /* num_channels channel steps in input order (ipstep_index.* fracts */
  double *step_fracts;  /* num_channels used to compute istep_index .* fracts */
  int *isteps_per_channel; /* num_channels */
  int *msteps_per_channel;
  int *istepped_xyz_ptrs;   /* num_channels */
  int *istep_index;        /* Current step number for each channel (for the sorted ordering) */
  int *ipstep_index;       /* step index in original order. */
  int *mm_istep_index;     /* maximum multiplie of istep_index that can 
			      occur in max_index power levels */
  int *icolors;            /* num_channels */
  int *imap_sort;          /* 2 * num_channels */
  int *refinement_intervals; /* max_scales */
  int *refinement_level;   /* num_channels, current refinement level for each channel */
  int *istep_interval;   /* Current interval (istep_index increment) for each channel */
  int *last_index;       /* num_channels, last istep index (used to be msteps for all */
  int *previous_index_in_bq; /* num_channels, if previous istep_index for a channel
				was not in the bounding quadrilateral, 1 if it was */
  int *unique_prime_factors; /* length is 4 * max_index if max_indx < 2310 
			                  5 * max_index if max_index < 30030*/
  int *upf_index;            /* length is max_index + 2 index of positions 
				of a numbers unique factor list.
				elements 0 and 1 are unused,
				elements i in [2:max_index]
				point to the start of integer i's unique prime factor
				list in unique_prime_factors vector.
				element [max_index+1] is one beyond max_index's list
				of prime factors. so Integer i's prime factors are
				in unique_prime_factors[upf_index[i]:upf_index[i+1]-1]
				for i in [2:max_index]. This array is used
			        to count the number of unique factors when 
			        building unique_prime_factors
			     */
  double vb_plane[4];
  double ub_plane[4];
  double u_target;
  double v_target;
  double test_rad;
  double u_a;
  double u_b;
  double v_a;
  double v_b;
  double mu_a;
  double mu_b;
  double mv_a;
  double mv_b;
  double recip_mu_a;
  double recip_mu_b;
  double recip_mv_a;
  double recip_mv_b;
  double xmin_c;
  double xmax_c;
  double ymin_c;
  double ymax_c;
  double xc_c;
  double yc_c;
  double xd_c;
  double yd_c;
  double x_top;
  double y_top;
  double z_top;
  double zmax;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double dvb;
  double dub;
  double y_dom_vb_min;
  double y_dom_vb_max;
  double y_dom_ub_min;
  double y_dom_ub_max;
  double x_dom_vb_min;
  double x_dom_vb_max;
  double x_dom_ub_min;
  double x_dom_ub_max;
  double z_dom_xsum;
  double z_dom_ysum;
  double z_dom_zsum;
  double y_dom_xsum;
  double y_dom_ysum;
  double y_dom_zsum;
  double x_dom_xsum;
  double x_dom_ysum;
  double x_dom_zsum;
  double va_delta; /* va_delta is the distance from 
		      (xmax_zmax,ymax_zmax) to (xc_zmax,yc_zmax)
		      it serves as an upper bound to the distance 
		      between line vb and line va */
  double ua_delta; /* ua_delta is the distance from
		      (xmax_zmax,ymax_zmax) to (xd_zmax,yd_zmax)
		      it serves as an upper bound to the distance 
		      between line ub and line ua */
  
  double min_scale;
  double three_mv_a;
  double three_mv_b;
  double min_fract;

  int num_channels;
  int nz;

  int ny;
  int nx;

  int max_index; /* New Maximum number of steps index = 1/min interval width */
  int number_scales;

  int task_msg_len;
  int task_msg_tag;

  FILE *lfp;
  FILE *ofp;
}
;
#endif
