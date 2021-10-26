#include "system_includes.h"
#include "state_struct.h"
#include "ts_sort.h"
#include "fill_isteps_per_channel.h"
#include "set_refinement_intervals.h"
#include "build_upf_list.h"
#include "plane_eq.h"
#include "init_loop_state.h"
int init_loop_state(struct state_struct *state,
	            int num_channels, double u_target, double v_target,
		    double test_rad) {
  /*
    Allocate and initialize the loop_state structure
    Called by: era_opt
    Calls:     ts_sort,fill_isteps_per_channel,set_refinement_intervals

    Arguments:     tmf       Description
    state          p*b       thous_opt state struct
                               lfp field is an input
			       use_multiscale field
			       loop_state field and subfields are output
    num_channels   i0i       Number of channels
    u_target       d0i       u-coordinate of the target circle.
    v_target       d0i       v-coordinate of the target circle.
    test_rad       d0i       radius of the target circle.

    isteps_per_channel
                   i1i       Number of steps per channel length is num_channels


    Uses the followingfields of state:
       lfp,
       number_scales,
       scales,
       test_radii_sq
       sbar_t_channels

    Sets the following fields of state
       loop_state and its subfields:
          refinement_intervals,
	  stepped_channel_xyz,
	  prev_xyz,
	  xyz_rps,
	  curr_xyz,
	  xyz_sort,
	  fracts,
	  pstep_fracts,
	  step_fracts,
	  istepped_xyz_ptrs,
	  icolors,
	  imap_sort,
	  msteps_per_channel,
	  refinement_level,
	  istep_interval,
	  last_index,
	  previous_index_in_bq,
	  num_channels,
	  nx,
	  ny,
	  nz,
	  u_a,
	  u_b,
	  v_a,
	  v_b,
	  mu_a,
	  mu_b,
	  mv_a,
	  mv_b
	  three_mv_a,
	  three_mv_b,
	  recip_mu_a,
	  recip_mu_b,
	  recip_mv_a,
	  recip_mv_b,
	  xmin_c,
	  ymin_c,
	  xmax_c,
	  ymax_c,
	  xc_c,
	  yc_c,
	  xd_c,
	  yd_c,
	  x_top,
	  y_top,
	  z_top,
	  zmax,
	  lfp
  */
  struct loop_state_struct *loop_state;
  struct loop_state_struct loop_state_i;
  double *ts_xyz;
  double *stepped_channel_xyz; /* 3 * sum(steps_per_channel) */
  double *prev_xyz; /* 3 * num_channels */
  double *xyz_rps; /* 3 * num_channels +3, reverse partials sums of 
		      tristimulus xyz */
  double *curr_xyz; /* 3 */
  double *xyz_sort; /* 6*num_channels) */
  double *fracts;   /* num_channels */
  double *pstep_fracts; /* num_channels */
  double *step_fracts; /* num_channels */
  double *sort_scratch;
  double *rps;
  double *xyzs;
  double *xyzstep;
  double *scales;  /* number_scales */
  double *test_radii_sq;
  double *vb_plane;
  double *ub_plane;
  double *tsv;
  double plane_points[9];
  double xmax_zmax;
  double ymax_zmax;
  double xmax_zmin;
  double ymax_zmin;
  double xd_zmax;
  double yd_zmax;
  double xc_zmax;
  double yc_zmax;
  double zmin;
  double plane_side;
  double dvb;
  double dub;
  double z_dom_xsum;
  double z_dom_ysum;
  double z_dom_zsum;
  double y_dom_xsum;
  double y_dom_ysum;
  double y_dom_zsum;
  double x_dom_xsum;
  double x_dom_ysum;
  double x_dom_zsum;
  double y_dom_vb_max;
  double y_dom_vb_min;
  double y_dom_ub_max;
  double y_dom_ub_min;
  double x_dom_vb_max;
  double x_dom_vb_min;
  double x_dom_ub_max;
  double x_dom_ub_min;

  double y_dom_vb_x;
  double y_dom_vb_y;
  double y_dom_vb_z;
  double y_dom_ub_x;
  double y_dom_ub_y;
  double y_dom_ub_z;
  double x_dom_vb_x;
  double x_dom_vb_y;
  double x_dom_vb_z;
  double x_dom_ub_x;
  double x_dom_ub_y;
  double x_dom_ub_z;
  double va_delta;
  double ua_delta;
  double new_va_delta;
  double new_ua_delta;

  int *isteps_per_channel;  /* num_channels */
  int *istepped_xyz_ptrs;   /* num_channels */
  int *istep_index;        /* Current step number for each channel in sorted order*/
  int *ipstep_index;       /* Current step number for each channel in raw input order */
  int *mm_istep_index;     /* maximum multiple of istep_index that can 
			      occur in max_index power levels */
  int *icolors;            /* num_channels */
  int *imap_sort;          /* mapping of sorted channel order to input
			     channel order 
 	                   */
  int *imap_sort2;
  int *msteps_per_channel;  /* num_channels */
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

  /* These only used if use_multiscale is 1 */
  int *refinement_intervals; /* max_scales */
  int *refinement_level;   /* num_channels, current refinement level for each channel */
  int *istep_interval;   /* Current interval (istep_index increment) for each channel */
  int *last_index;       /* num_channels, last istep index (used to be msteps for all */
  int *previous_index_in_bq; /* num_channels, if previous istep_index for a channel
				was not in the bounding quadrilateral, 1 if it was */

  int64_t ask_for;
  int64_t one_l;
  int64_t usage;
  int64_t sizeof_double;
  int64_t sizeof_int;

  double min_scale;
  double u_a;
  double u_b;
  double v_a;
  double v_b;
  double mu_a;
  double mu_b;
  double mv_a;
  double mv_b;
  double xmin_c;
  double xmax_c;
  double ymin_c;
  double ymax_c;
  double x_top;
  double y_top;
  double z_top;
  double zmax;
  double zmax_2;
  double x;
  double y;
  double z;
  double fract;
  double fract_j;
  double xc_c;
  double yc_c;
  double xd_c;
  double yd_c;
  double min_fract;

  int nz;
  int ny;

  int nx;
  int success;

  int isteps_total;
  int i;

  int j;
  int k;

  int ipos[3];
  int ig;

  int ic;
  int ifield;

  int i3;
  int number_scales;

  int use_multiscale;
  int max_index;

  int first_interval;
  int extra_log_info;

  int colinear;
  int i3nz;

  int i3nzpny;
  int mpf;

  int use_bruteforce;
  int ipad;

  FILE *lfp;
  FILE *efp;

  success            = 1;
  usage              = (int64_t)0;
  lfp                = state->lfp;
  scales             = state->scales;
  test_radii_sq      = state->test_radii_sq;
  use_multiscale     = state->use_multiscale;
  use_bruteforce     = state->use_bruteforce;
  number_scales      = state->number_scales;
  scales             = state->scales;
  isteps_per_channel = state->isteps_per_channel;
  extra_log_info     = state->extra_log_info;
  sizeof_double      = (int64_t)(sizeof(double));
  sizeof_int         = (int64_t)(sizeof(int));
  one_l              = (int64_t)1;
  refinement_intervals = NULL;
  stepped_channel_xyz  = NULL;
  prev_xyz             = NULL;
  xyz_rps              = NULL;
  curr_xyz             = NULL;
  xyz_sort             = NULL;
  fracts               = NULL;
  pstep_fracts         = NULL;
  step_fracts          = NULL;
  istepped_xyz_ptrs    = NULL;
  icolors              = NULL;
  imap_sort            = NULL;
  msteps_per_channel   = NULL;
  refinement_level     = NULL;
  istep_interval       = NULL;
  last_index           = NULL;
  previous_index_in_bq = NULL;
  unique_prime_factors = NULL;
  upf_index            = NULL;

  ask_for    = (int64_t)sizeof(loop_state_i);
  usage      += ask_for;
  loop_state = (struct loop_state_struct *)calloc(one_l,ask_for);
  if (loop_state == NULL) {
    success = 0;
    if (lfp) {
      fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for loop_state struct \n",ask_for);
      fflush(lfp);
    }
  }
  isteps_total = 0;
  if (success) {
    state->loop_state              = loop_state;
    loop_state->num_channels       = num_channels;
    loop_state->u_target           = u_target;
    loop_state->v_target           = v_target;
    loop_state->test_rad           = test_rad;
    loop_state->number_scales      = number_scales;
    loop_state->isteps_per_channel = isteps_per_channel;
    vb_plane                       = (double*)&loop_state->vb_plane[0];
    ub_plane                       = (double*)&loop_state->ub_plane[0];
  }
  isteps_total = 0;
  if (success) {
    if (use_multiscale != 1) {
      fill_isteps_per_channel(state);
      isteps_total = 0;
      max_index = 0;
      for (i=0;i<num_channels;i++) {
	if (isteps_per_channel[i] > max_index) {
	  max_index = isteps_per_channel[i];
	}
	isteps_total += isteps_per_channel[i];
      }
      max_index = max_index - 1;
    } else {
      ask_for              = (int64_t)(number_scales * sizeof_int);
      usage                += ask_for;
      refinement_intervals = (int*)calloc(one_l,ask_for);
      if (refinement_intervals == NULL) {
	success = 0;
	if (lfp) {
	  fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for"
		   " refinement_intervals.\n",ask_for);
	  fflush(lfp);
	}
      }
      if (success) {
	loop_state->refinement_intervals = refinement_intervals;
	set_refinement_intervals(number_scales,num_channels,scales,
				 &min_scale,refinement_intervals,
				 &max_index,&isteps_total,lfp);
				 
	loop_state->max_index = max_index;
	loop_state->min_scale = min_scale;
      }
    }
    if (isteps_total == 0) {
      success = 0;
    }
  }
  /*
    Now that we have max_index we can allocate unique_prime_factors and
    upf_index vectors. note that 2*3*5*7 = 210 and 2*3*5*7*11 = 2310,
    so for max_index < 2310 numbers will have at most 4 unique prime factors.
    For numbers < 2310*13 = 30030, numbers will have at most 5 unqiue prime factors.
  */
  if (success) {
    if (max_index < 2310) {
      mpf = 4;
    } else {
      mpf = 5;
      if (max_index >= 30030) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"init_loop_state: Code not set to handle duplicate ratio"
                      "removal for max_index >= 30030. "
		  "\nContact code Maintainer Doug Baxter\n");
	  fflush(lfp);
	}
      }
    }
  }
  if (success) {
    ask_for = (int64_t)(max_index * mpf * sizeof(int));
    usage += ask_for;
    unique_prime_factors = (int*)calloc(one_l,ask_for);
    if (unique_prime_factors == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"init_loop_state: Could not allcoate %ld bytes for "
		"unique_prime_factors\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    ask_for = (int64_t)((max_index + 2) * sizeof(int));
    usage   += ask_for;
    upf_index = (int *)calloc(one_l,ask_for);
    if (upf_index == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"init_loop_state: Could not allcoate %ld bytes for "
		"upf_index\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    ask_for = (int64_t)(3*isteps_total * sizeof(double));
    usage += ask_for;
    stepped_channel_xyz = (double *)calloc(one_l,ask_for);
    if (stepped_channel_xyz == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for stepped_channel_xyz \n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->stepped_channel_xyz = stepped_channel_xyz;
    ask_for = (int64_t)(3*num_channels * sizeof_double);
    usage += ask_for;
    prev_xyz = (double *)calloc(one_l,ask_for);
    if (prev_xyz == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for prev_xyz \n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->prev_xyz = prev_xyz;
    ask_for = (int64_t)(3*(num_channels+1) * sizeof_double);
    usage += ask_for;
    xyz_rps = (double *)calloc(one_l,ask_for);
    if (xyz_rps == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for xyz_rps \n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->xyz_rps = xyz_rps;
    ask_for = (int64_t)(4*sizeof_double);
    usage += ask_for;
    curr_xyz = (double *)calloc(one_l,ask_for);
    if (curr_xyz == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for curr_xyz\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->curr_xyz = curr_xyz;
    ask_for = (int64_t)(6*num_channels*sizeof_double);
    usage += ask_for;
    xyz_sort = (double *) calloc(one_l,ask_for);
    if (xyz_sort == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for xyz_sort\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->xyz_sort = xyz_sort;
    ask_for = (int64_t)(num_channels * sizeof_double);
    usage += ask_for;
    fracts = (double *) calloc(one_l,ask_for);
    if (fracts == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for fracts\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->fracts = fracts;
    ask_for = (int64_t)(num_channels * sizeof_double);
    usage += ask_for;
    pstep_fracts = (double *) calloc(one_l,ask_for);
    if (pstep_fracts == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for pstep_fracts\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->pstep_fracts = pstep_fracts;
    ask_for = (int64_t)(num_channels * sizeof_double);
    usage += ask_for;
    step_fracts = (double *) calloc(one_l,ask_for);
    if (step_fracts == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for step_fracts\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->step_fracts = step_fracts;
    loop_state->isteps_per_channel = isteps_per_channel;
    ask_for = (int64_t) (num_channels * sizeof(int));
    usage += ask_for;
    istepped_xyz_ptrs = (int *) calloc(one_l,ask_for);
    if (istepped_xyz_ptrs == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld "
		 "bytes for istepped_xyz_ptrs\n",ask_for);
	fflush(lfp);
      }
    } else {
      loop_state->istepped_xyz_ptrs = istepped_xyz_ptrs;
    }
  }
  if (success) {
    ask_for = (int64_t) (num_channels * sizeof(int));
    usage += ask_for;
    istep_index = (int *)calloc(one_l,ask_for);
    if (istep_index == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld "
		 "bytes for istep_index\n",ask_for);
	fflush(lfp);
      }
    } else {
      loop_state->istep_index = istep_index;
    }
  }
  if (success) {
    ask_for = (int64_t) (num_channels * sizeof(int));
    usage += ask_for;
    ipstep_index = (int *)calloc(one_l,ask_for);
    if (ipstep_index == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld "
		 "bytes for ipstep_index\n",ask_for);
	fflush(lfp);
      }
    } else {
      loop_state->ipstep_index = ipstep_index;
    }
  }
  if (success) {
    ask_for = (int64_t) (num_channels * sizeof(int));
    usage += ask_for;
    mm_istep_index = (int *)calloc(one_l,ask_for);
    if (mm_istep_index == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for "
		 "mm_istep_index\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->mm_istep_index = mm_istep_index;
    ask_for = (int64_t) (num_channels * sizeof(int));
    usage += ask_for;
    icolors = (int *)calloc(one_l,ask_for);
    if (icolors == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for icolors\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->icolors = icolors;
    ask_for = (int64_t) (2 * num_channels * sizeof(int));
    usage += ask_for;
    imap_sort = (int *)calloc(one_l,ask_for);
    if (imap_sort == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for imap_sort\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->imap_sort = imap_sort;
    ask_for = (int64_t)(num_channels * sizeof(int));
    usage += ask_for;
    msteps_per_channel = (int *)calloc(one_l,ask_for);
    if (msteps_per_channel == NULL) {
      success = 0;
      if (lfp) {
	fprintf (lfp,"init_loop_state: Error could not allocate %ld bytes for msteps_per_channel\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    loop_state->msteps_per_channel = msteps_per_channel;
    /*
      Check validity of isteps_per_channel input.
    */
    for (i=0;i<num_channels;i++) {
      if (isteps_per_channel[i] < 2) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"init_loop_state: isteps_per_channel must be at least 2."
		    " isteps_per_channel[%d] = %d\n",i,isteps_per_channel[i]);
	  fflush(lfp);
	}
      }
    }
    /*
      Allocate arrays for multiscale combination processing.
    */
    ask_for  = (int64_t)(num_channels * sizeof_int);
    usage   += ask_for;
    refinement_level = (int*) calloc(one_l,ask_for);
    if (refinement_level == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"init_loop_state: Error could not allocate %ld bytes for "
		"refinement_level\n",ask_for);
	fflush(lfp);
      }
    }
    if (success) {
      loop_state->refinement_level = refinement_level;
      ask_for  = (int64_t)(num_channels * sizeof_int);
      usage   += ask_for;
      istep_interval = (int*) calloc(one_l,ask_for);
      if (istep_interval == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"init_loop_state: Error could not allocate %ld bytes for "
		  "istep_interval\n",ask_for);
	  fflush(lfp);
	}
      }
    }
    if (success) {
      loop_state->istep_interval = istep_interval;
      ask_for  = (int64_t)(num_channels * sizeof_int);
      usage   += ask_for;
      last_index = (int*) calloc(one_l,ask_for);
      if (last_index == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"init_loop_state: Error could not allocate %ld bytes for "
		  "last_index\n",ask_for);
	  fflush(lfp);
	}
      }
    }
    if (success) {
      loop_state->last_index = last_index;
      ask_for  = (int64_t)(num_channels * sizeof_int);
      usage   += ask_for;
      previous_index_in_bq = (int*) calloc(one_l,ask_for);
      if (previous_index_in_bq == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"init_loop_state: Error could not allocate %ld bytes for "
		  "previous_index_in_bq\n",ask_for);
	  fflush(lfp);
	}
      }
    }
    if (success) {
      loop_state->previous_index_in_bq = previous_index_in_bq;
    }
  }
  if (success) {
    ts_xyz = state->sbar_t_channels;
    if (extra_log_info) {
      if (lfp) {
	fprintf(lfp,"Tristimulus xyz values\n");
	i3 = 0;
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%le\t%le\t%le\n",ts_xyz[i3],ts_xyz[i3+1],ts_xyz[i3+2]);
	  i3 += 3;
	}
	fflush(lfp);
      }
    }
    /*
      Count the channel colors (nx,ny,nz)
    */
    nx = 0;
    ny = 0;
    nz = 0;
 
    j = 0;
    for (i=0;i<num_channels;i++) {
      x = ts_xyz[j];
      y = ts_xyz[j+1];
      z = ts_xyz[j+2];
      if ((z > x) && (z > y)) {
	icolors[i] = 0;
	nz = nz + 1;
      } else if (x > y) {
	icolors[i] = 2;
	nx = nx + 1;
      } else {
	icolors[i] = 1;
	ny = ny + 1;
      }
      j = j + 3;
    }
    loop_state->nx = nx;
    loop_state->ny = ny;
    loop_state->nz = nz;
    /*
      need to check for least one channel
      with dominant tsv coordinate for 
      each coordinate.
    */
    if (lfp) {
      fprintf(lfp,"nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
      fflush(lfp);
    }
    if ((nx == 0) || (ny == 0) || (nz == 0)) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"init_loop_state: Error, you must have at least"
		" one channel with a dominant TSV coordinate for each"
		" coordinate, nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      Group points of 1 color together
      And put the z points first followed by y then x.
    */
    ipos[0] = 0;
    ipos[1] = nz;
    ipos[2] = nz + ny;
    j = 0;
    for (i=0;i<num_channels;i++) {
      ic = icolors[i];
      ig = ipos[ic];
      imap_sort[ig] = i;
      ig = ig * 3;
      xyz_sort[ig]   = ts_xyz[j];
      xyz_sort[ig+1] = ts_xyz[j+1];
      xyz_sort[ig+2] = ts_xyz[j+2];
      j = j + 3;
      ipos[ic] += 1;
    }
    ipos[0] = 0;
    ipos[1] = nz;
    ipos[2] = nz + ny;
    if (extra_log_info) {
      if (lfp) {
	fprintf(lfp,"grouped Tristimulus xyz values\n");
	fprintf(lfp,"nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
	i3 = 0;
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%le\t%le\t%le\n",xyz_sort[i3],xyz_sort[i3+1],
		  xyz_sort[i3+2]);
	  i3 += 3;
	}
	fprintf (lfp,"\nimap_sort = ");
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%d\t",imap_sort[i]);
	}
	fprintf(lfp,"\n");
	fflush(lfp);
      }
    }
    /*
      Now sort each color in descending order on its coordinate strength.
    */
    sort_scratch = &xyz_sort[3*num_channels];
    imap_sort2   = &imap_sort[num_channels];

    ifield = 2;
    ts_sort(nz,ifield,xyz_sort,imap_sort,sort_scratch,imap_sort2);

    if (extra_log_info) {
      if (lfp) {
	fprintf(lfp,"Tristimulus sorted z xyz values\n");
	fprintf(lfp,"nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
	i3 = 0;
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%le\t%le\t%le\n",xyz_sort[i3],xyz_sort[i3+1],
		  xyz_sort[i3+2]);
	  i3 += 3;
	}
	fprintf (lfp,"\nimap_sort = ");
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%d\t",imap_sort[i]);
	}
	fprintf(lfp,"\n");
	fflush(lfp);
      }
    }
    ic = ipos[1];
    ig = ic * 3;
    ifield = 1;
    ts_sort(ny,ifield,&xyz_sort[ig],&imap_sort[ic],sort_scratch,imap_sort2);
    if (extra_log_info) {
      if (lfp) {
	fprintf(lfp,"Tristimulus sorted y xyz values\n");
	fprintf(lfp,"nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
	i3 = 0;
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%le\t%le\t%le\n",xyz_sort[i3],xyz_sort[i3+1],xyz_sort[i3+2]);
	  i3 += 3;
	}
	fprintf (lfp,"\nimap_sort = ");
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%d\t",imap_sort[i]);
	}
	fprintf(lfp,"\n");
	fflush(lfp);
      }
    }
    ic = ipos[2];
    ig = ic * 3;
    ifield = 0;
    ts_sort(nx,ifield,&xyz_sort[ig],&imap_sort[ic],sort_scratch,imap_sort2);
    if (extra_log_info) {
      if (lfp) {
	fprintf(lfp,"sorted Tristimulus xyz values\n");
	fprintf(lfp,"nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
	i3 = 0;
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%le\t%le\t%le\n",xyz_sort[i3],xyz_sort[i3+1],
		  xyz_sort[i3+2]);
	  i3 += 3;
	}
	fprintf (lfp,"\nimap_sort = ");
	for (i=0;i<num_channels;i++) {
	  fprintf(lfp,"%d\t",imap_sort[i]);
	}
	fprintf(lfp,"\n");
	fflush(lfp);
      }
    }
    /*
      Now it is helpful in determining the number of bins for
      sorting the colors to have the xmax,ymax,zmax values within
      each set of channels.
    */
    z_dom_xsum = 0.0;
    z_dom_ysum = 0.0;
    z_dom_zsum = 0.0;
    tsv       = &xyz_sort[0];
    for (i=0;i<nz;i++) {
      z_dom_xsum += tsv[0];
      z_dom_ysum += tsv[1];
      z_dom_zsum += tsv[2];
      tsv += 3;
    }
    i3nz       = 3*nz;
    y_dom_xsum = 0.0;
    y_dom_ysum = 0.0;
    y_dom_zsum = 0.0;
    tsv       = &xyz_sort[i3nz];
    for (i=nz;i<nz+ny;i++) {
      y_dom_xsum += tsv[0];
      y_dom_ysum += tsv[1];
      y_dom_zsum += tsv[2];
      tsv += 3;
    }
    i3nzpny  = 3*(nz+ny);
    x_dom_xsum = 0.0;
    x_dom_ysum = 0.0;
    x_dom_zsum = 0.0;
    tsv       = &xyz_sort[i3nzpny];
    for (i=ny+nz;i<nz+ny+nx;i++) {
      x_dom_xsum += tsv[0];
      x_dom_ysum += tsv[1];
      x_dom_zsum += tsv[2];
      tsv += 3;
    }
    /*
      Compute wedge slopes
      u_target = .2357;
      v_target = .5113;

      test_rad = 0.005;
    */

    u_a = u_target - test_rad;
    u_b = u_target + test_rad;

    v_a = v_target - test_rad;
    v_b = v_target + test_rad;

    mu_a = (4.0-u_a)/(15.0 * u_a);
    mu_b = (4.0-u_b)/(15.0 * u_b);
    mv_a = v_a/(9.0-(15.0*v_a));
    mv_b = v_b/(9.0-(15.0*v_b));
/*

    (xmin,ymin) is at the intersection of  the u_a and v_a lines
  xmin = ((3mv_a + .2)/(mu_a - mv_a)) * z  = xmin_c * z
  ymin = mv_a * xmin + 3mv_a * z  =          ymin_c * z

       = (mv_a*((3mv_a + .2)/(mu_a-mv_a)) + 3mv_a)*z

    (xmax,ymax) is at the intersection of the u_b and v_b lines
  xmax = ((3mv_b + .2)/(mu_b -mv_b)) * z   = xmax_c * z
  ymax = mv_b * xmax + 3mv_b * z           = ymax_c * z

    (xc,yc) will be the intersection of the u_b and v_a lines
    xc = ((3mv_a + .2)/(mu_b - mv_a))z  = xc_c * z
    yc = ((3mv_b + .2)*mv_a/(mu_b-mv_a))z = yc_c & z

    (xc,yd) will be the intersection of the u_a and v_b lines
    xd = (.2 + 3mv_b)/(mu_a - mv_b)

*/
    xmin_c = ((3.0*mv_a) + 0.2)/(mu_a - mv_a);
    ymin_c = mv_a * (xmin_c + 3.0);
    xmax_c = ((3.0*mv_b) + 0.2)/(mu_b - mv_b);
    ymax_c = mv_b * (xmax_c + 3.0);
  /*
  leeway = 0.01;
  mratio = 1.0-leeway;
  pratio = 1.0+leeway;
  xmin_c = xmin_c * mratio;
  ymin_c = ymin_c * mratio;
  xmax_c = xmax_c * pratio;
  ymax_c = ymax_c * pratio;
  */
    xc_c   = ((3.0*mv_a) + 0.2)/(mu_b - mv_a);
    yc_c   = mv_a * (xc_c + 3.0);
    xd_c   = ((3.0*mv_b) + 0.2)/(mu_a - mv_b);
    yd_c   = mv_b * (xd_c + 3.0);

    rps = &xyz_rps[(num_channels)*3];
    xyzs = &xyz_sort[(num_channels-1)*3];
    rps[0] = 0.0;
    rps[1] = 0.0;
    rps[2] = 0.0;
    rps = rps - 3;    /* Caution address arithmetic */
    for (i= num_channels - 1;i>=0;i--) {
      rps[0] = rps[3] + xyzs[0];
      rps[1] = rps[4] + xyzs[1];
      rps[2] = rps[5] + xyzs[2];
      rps = rps - 3;    /* Caution address arithmetic */
      xyzs = xyzs - 3;  /* Caution address arithmetic */
    }
    rps = rps + 3;      /* Caution address arithmetic */
    x_top = rps[0];
    y_top = rps[1];
    z_top = rps[2];
    /*
      z <=  x_tops *(mu_a - mv_a) /(3mv_a + .2)

      and z < = y_tops/(mv_a * ((3mv_a + .2)/(mu_a-mv_a)) + 3mv_a)
    */
    zmax  =  x_top / xmin_c;
    zmax_2  =  y_top / ymin_c;
    if (zmax_2 < zmax) {
      zmax = zmax_2;
    }
    /*
      Now we want to compute the coordinates for the
      vb plane and the ub plane. To do this we need three points on each plane 
      For the vb plane we will use xmax,ymax,zmax, xd,yd,zmax, and 
      xmax,ymax,zmin
      For the ub plane we will use xmax,ymax,zmax, xc,yc,zmax, xmax,ymax,zmin

      We also want dvb the maximum vertical distance between the
      v_b and v_a lines and dub the maximumu horizontal distance between the
      u_b and u_a lines. Both of theses work out to be about .003
    */

    min_fract = 1.0/((double)max_index);
    i3        = (nz-1)*3;
    zmin      = xyz_sort[i3+2]*min_fract;

    xmax_zmax = xmax_c * zmax;
    ymax_zmax = ymax_c * zmax;
    xmax_zmin = xmax_c * zmin;
    ymax_zmin = ymax_c * zmin;
    xd_zmax   = xd_c * zmax;
    yd_zmax   = yd_c * zmax;
    xc_zmax   = xc_c * zmax;
    yc_zmax   = yc_c * zmax;
    va_delta = (xmax_zmax - xc_zmax) * (xmax_zmax - xc_zmax);
    va_delta = va_delta + ((ymax_zmax - yc_zmax) * (ymax_zmax - yc_zmax));
    va_delta = sqrt(va_delta);
    
    ua_delta = (xmax_zmax - xd_zmax) * (xmax_zmax - xd_zmax);
    ua_delta = ua_delta + ((ymax_zmax - yd_zmax) * (ymax_zmax - yd_zmax));
    ua_delta = sqrt(ua_delta);

    dvb = (mv_b-mv_a)*(xmax_zmax+(3*zmax));
    dub = (ymax_zmax+(.2*zmax))*(mu_a-mu_b)/(mu_a*mu_b);

    plane_points[0] = xmax_zmin;
    plane_points[1] = ymax_zmin;
    plane_points[2] = zmin;
    plane_points[3] = xmax_zmax;
    plane_points[4] = ymax_zmax;
    plane_points[5] = zmax;
    plane_points[6] = xd_zmax;
    plane_points[7] = yd_zmax;
    plane_points[8] = zmax;

    plane_eq(plane_points,vb_plane,&colinear);
    if (colinear) {
      if (lfp) {
	success = 0;
	fprintf(lfp,"init_loop_state: vb_plane points were colinear\n");
	fflush(lfp);
      }
    } else {
      /*
	Strangely enough the normalization always get the signs reversed
	for vb_plane based on which side of the plane xc is (is should be on the negative side).
	So we just change the sign of the vb_plane elements.
      */
      for (i=0;i<3;i++) {
	vb_plane[i] = -vb_plane[i];
      }
      plane_points[0] = xmax_zmax;
      plane_points[1] = ymax_zmax;
      plane_points[2] = zmax;
      plane_points[3] = xmax_zmin;
      plane_points[4] = ymax_zmin;
      plane_points[5] = zmin;
      plane_points[6] = xc_zmax;
      plane_points[7] = yc_zmax;
      plane_points[8] = zmax;
      plane_eq(plane_points,ub_plane,&colinear);
      if (colinear) {
	if (lfp) {
	  success = 0;
	  fprintf(lfp,"init_loop_state: ub_plane points were colinear\n");
	  fflush(lfp);
	}
      }
    }
  }
  if (success) {
    y_dom_vb_max = 0.0;
    y_dom_vb_min = 0.0;
    y_dom_vb_x = vb_plane[0] * y_dom_xsum;
    if (y_dom_vb_x > 0.0) {
      y_dom_vb_max = y_dom_vb_x;
    } else {
      y_dom_vb_min = y_dom_vb_x;
    }
    y_dom_vb_y = vb_plane[1] * y_dom_ysum;
    if (y_dom_vb_y > 0.0) {
      y_dom_vb_max += y_dom_vb_y;
    } else {
      y_dom_vb_min += y_dom_vb_y;
    }
    y_dom_vb_z = vb_plane[2] * y_dom_zsum;
    if (y_dom_vb_z > 0.0) {
      y_dom_vb_max += y_dom_vb_z;
    } else {
      y_dom_vb_min += y_dom_vb_z;
    }
    y_dom_ub_max = 0.0;
    y_dom_ub_min = 0.0;
    y_dom_ub_x = ub_plane[0] * y_dom_xsum;
    if (y_dom_ub_x > 0.0) {
      y_dom_ub_max = y_dom_ub_x;
    } else {
      y_dom_ub_min = y_dom_ub_x;
    }
    y_dom_ub_y = ub_plane[1] * y_dom_ysum;
    if (y_dom_ub_y > 0.0) {
      y_dom_ub_max += y_dom_ub_y;
    } else {
      y_dom_ub_min += y_dom_ub_y;
    }
    y_dom_ub_z = ub_plane[2] * y_dom_zsum;
    if (y_dom_ub_z > 0.0) {
      y_dom_ub_max += y_dom_ub_z;
    } else {
      y_dom_ub_min += y_dom_ub_z;
    }

    x_dom_vb_max = 0.0;
    x_dom_vb_min = 0.0;
    x_dom_vb_x = vb_plane[0] * x_dom_xsum;
    if (x_dom_vb_x > 0.0) {
      x_dom_vb_max = x_dom_vb_x;
    } else {
      x_dom_vb_min = x_dom_vb_x;
    }
    x_dom_vb_y = vb_plane[1] * x_dom_ysum;
    if (x_dom_vb_y > 0.0) {
      x_dom_vb_max += x_dom_vb_y;
    } else {
      x_dom_vb_min += x_dom_vb_y;
    }
    x_dom_vb_z = vb_plane[2] * x_dom_zsum;
    if (x_dom_vb_z > 0.0) {
      x_dom_vb_max += x_dom_vb_z;
    } else {
      x_dom_vb_min += x_dom_vb_z;
    }
    x_dom_ub_max = 0.0;
    x_dom_ub_min = 0.0;
    x_dom_ub_x = ub_plane[0] * x_dom_xsum;
    if (x_dom_ub_x > 0.0) {
      x_dom_ub_max = x_dom_ub_x;
    } else {
      x_dom_ub_min = x_dom_ub_x;
    }
    x_dom_ub_y = ub_plane[1] * x_dom_ysum;
    if (x_dom_ub_y > 0.0) {
      x_dom_ub_max += x_dom_ub_y;
    } else {
      x_dom_ub_min += x_dom_ub_y;
    }
    x_dom_ub_z = ub_plane[2] * x_dom_zsum;
    if (x_dom_ub_z > 0.0) {
      x_dom_ub_max += x_dom_ub_z;
    } else {
      x_dom_ub_min += x_dom_ub_z;
    }
  }

  if (success) {
    build_upf_list(max_index,mpf,upf_index,unique_prime_factors);
    if (extra_log_info) {
      if (lfp) {
	fprintf(lfp,"test_rad = %le\n",test_rad);
	fprintf(lfp,"u_a = %le, u_target = %le, u_b = %le\n",
		u_a,u_target,u_b);
	fprintf(lfp,"v_a = %le, v_target = %le, v_b = %le\n",
		v_a,v_target,v_b);
	fprintf(lfp,"mu_a = %le, mu_b = %le\n",mu_a,mu_b);
	fprintf(lfp,"mv_a = %le, mv_b = %le\n",mv_a,mv_b);
	fprintf(lfp,"xmin_c = %le, xmax_c = %le\n",xmin_c,xmax_c);
	fprintf(lfp,"ymin_c = %le, ymax_c = %le\n",ymin_c,ymax_c);
	fprintf(lfp,"xc_c = %le, yc_c = %le\n",xc_c,yc_c);
	fprintf(lfp,"xd_c = %le, yd_c = %le\n",xd_c,yd_c);
	fprintf(lfp,"x_top = %le, y_top = %le, z_top = %le\n",
		x_top,y_top,z_top);
	fprintf(lfp,"zmax = %le\n",zmax);
	fprintf(lfp,"zmin = %le\n",zmin);
	fprintf(lfp,"vb_plane = %le, %le, %le, %le\n",vb_plane[0],
		vb_plane[1],vb_plane[2],vb_plane[3]);
	fprintf(lfp,"ub_plane = %le, %le, %le, %le\n",ub_plane[0],
		ub_plane[1],ub_plane[2],ub_plane[3]);
	plane_side    = (xc_zmax*vb_plane[0]) + (yc_zmax*vb_plane[1]) +
	                (zmax*vb_plane[2]) +  vb_plane[3];
	new_va_delta = -plane_side;
	/*
	  NB The absolute value of the plane side above should be va_delta.
	  Need to try that out. Should make for a narrower strip of candidates.
	*/
	fprintf(lfp,"plane side of xc_zmax,yc_zmax,zmax,vb_plane = %le\n",
		plane_side);
	fprintf(lfp,"old va_delta = %le, new va_delta = %le",va_delta,new_va_delta);
	va_delta = new_va_delta;
	plane_side    = (xmax_zmax*vb_plane[0]) + 
	                ((ymax_zmax+1.0)*vb_plane[1]) +
	                (zmax*vb_plane[2]) +  vb_plane[3];
	fprintf(lfp,"plane side of xmax_zmax,(ymax_zmax+1),zmax,vb_plane = %le\n",
		plane_side);
	plane_side    = (xd_zmax*ub_plane[0]) + (yd_zmax*ub_plane[1]) +
	                (zmax*ub_plane[2]) +  ub_plane[3];
	new_ua_delta = plane_side;
	fprintf(lfp,"old ua_delta = %le, new va_delta = %le",ua_delta,new_ua_delta);
	ua_delta = new_ua_delta;
	/*
	  ua_delta = plane_side;
	  The absolute value of the plane_size above should probablby be ua_delta.
	  Need to try that out. Should make for a narrower strip of candidates.
	*/
	fprintf(lfp,"plane side of xd_zmax,yd_zmax,zmax,ub_plane = %le\n",
		plane_side);
	plane_side    = (xmax_zmax*ub_plane[0]) + 
	                ((ymax_zmax-1.0)*ub_plane[1]) +
	                (zmax*ub_plane[2]) +  ub_plane[3];
	fprintf(lfp,"plane side of xmax_zmax,(ymax_zmax-1),zmax,ub_plane = %le\n",
		plane_side);
	fprintf(lfp,"dvb = %le,  dub = %le\n",dvb,dub);

	fprintf(lfp,"z_dom_xsum  = %le, z_dom_ysum  = %le, z_dom_zsum  = %le\n",
		z_dom_xsum,z_dom_ysum,z_dom_zsum);
	fprintf(lfp,"y_dom_xsum = %le, y_dom_ysum = %le, y_dom_zsum = %le\n",
		y_dom_xsum,y_dom_ysum,y_dom_zsum);
	fprintf(lfp,"x_dom_xsum   = %le, x_dom_ysum   = %le, x_dom_zsum   = %le\n",
		x_dom_xsum,x_dom_ysum,x_dom_zsum);

	fprintf(lfp,"y_dom_vb_max = %le\n",y_dom_vb_max);
	fprintf(lfp,"y_dom_vb_min = %le\n",y_dom_vb_min);
	fprintf(lfp,"y_dom_ub_max = %le\n",y_dom_ub_max);
	fprintf(lfp,"y_dom_ub_min = %le\n",y_dom_ub_min);
	fprintf(lfp,"x_dom_vb_max   = %le\n",x_dom_vb_max);
	fprintf(lfp,"x_dom_vb_min   = %le\n",x_dom_vb_min);
	fprintf(lfp,"x_dom_ub_max   = %le\n",x_dom_ub_max);
	fprintf(lfp,"x_dom_ub_min   = %le\n",x_dom_ub_min);
	fflush(lfp);
      }
    }
    /*
      Need to permute the isteps_per_channel to the new ordering.
    */
    if (use_multiscale != 1) {
      /*
	Need to permute the isteps_per_channel to the new ordering.
      */
      for (i = 0;i<num_channels;i++) {
	msteps_per_channel[i] = isteps_per_channel[imap_sort[i]];
      }
    } else {
      for (i = 0;i<num_channels;i++) {
	msteps_per_channel[i] = max_index + 1;
      }
    }
    if (extra_log_info) {
      if (lfp) {
	fprintf(lfp,"use_multiscale = %d\n",use_multiscale);
	fprintf(lfp,"msteps_per_channel = ");
	for (i = 0;i<num_channels;i++) {
	  fprintf(lfp,"%5d",msteps_per_channel[i]);
	}
	fprintf(lfp,"\n");
	if (use_multiscale == 1) {
	  fprintf(lfp,"number_scales = %d\n",number_scales);
	  fprintf(lfp,"min_scale     = %le\n",min_scale);
	  fprintf(lfp,"max_index     = %d\n",max_index);
	  fprintf(lfp,"refinement_intervals = ");
	  for (i=0;i<number_scales;i++) {
	    fprintf(lfp,"%5d",refinement_intervals[i]);
	  }
	  fprintf(lfp,"\n");
	  fprintf(lfp,"scales = ");
	  for (i=0;i<number_scales;i++) {
	    fprintf(lfp," %le",scales[i]);
	  }
	  fprintf(lfp,"\n");
	  fprintf(lfp,"\nistep_index[0:%d] , Rf, Rg, dist, nzc\n",num_channels);
	  fflush(lfp);
	}
	fprintf(lfp,"test_radii_sq = ");
	for (i=1;i<=num_channels;i++) {
	  fprintf(lfp,"%le ",test_radii_sq[i]);
	}
	fprintf(lfp,"\n");
	fflush(lfp);
      }
    }
    /*
      Now need to compute the pointers to the stepped values of each channel
      by forming the partial sums of 3 * msteps_per_channel.
    */
    istepped_xyz_ptrs[0] = 0;
    for (i=0;i<num_channels-1;i++) {
      istepped_xyz_ptrs[i+1] = istepped_xyz_ptrs[i] + (msteps_per_channel[i]*3);
    }
    /*
      Now we compute and store the stepped values of each channel.
    */
    xyzstep = stepped_channel_xyz;
    xyzs    = xyz_sort;
    for (i=0;i<num_channels;i++) {
      k = msteps_per_channel[i] - 1;
      fract = 1.0/((double)k);
      fracts[i] = fract;
      /*
	Do the zero step explicitly.
      */
      xyzstep[0] = 0.0;
      xyzstep[1] = 0.0;
      xyzstep[2] = 0.0;
      xyzstep += 3; /* Caution address arithmetic */
      x = xyzs[0];
      y = xyzs[1];
      z = xyzs[2];
      xyzs += 3; /* Caution address arithmetic */
      for (j=1;j<k;j++) {
	fract_j = fract * ((double)j);
	xyzstep[0] = fract_j*x;	
	xyzstep[1] = fract_j*y;
	xyzstep[2] = fract_j*z;
	xyzstep += 3; /* Caution address arithmetic */
      }
      /*
	DO the 100% step explicitly
      */
      xyzstep[0] = x;
      xyzstep[1] = y;
      xyzstep[2] = z;
      xyzstep += 3; /* Caution address arithmetic */
    }
    /*
      Zero out the prev_xyz array.
    */
    for (i=0;i<3*num_channels;i++) {
      prev_xyz[i] = 0.0;
    }
    /*
      Zero out the istep_index array.
    */
    for (i=0;i<num_channels;i++) {
      istep_index[i] = 0;
    }
    if (use_multiscale == 1) {
      for (i=0;i<num_channels;i++) {
	last_index[i] = max_index;
      }
      for (i=0;i<num_channels;i++) {
	refinement_level[i] = 0;
      }
      for (i=0;i<num_channels;i++) {
	previous_index_in_bq[i] = 0;
      }
      first_interval = refinement_intervals[0];
      for (i=0;i<num_channels;i++) {
	istep_interval[i] = first_interval;
      }
    } else {
      for (i=0;i<num_channels;i++) {
	last_index[i] = msteps_per_channel[i];
      }
      for (i=0;i<num_channels;i++) {
	refinement_level[i] = 0;
      }
      for (i=0;i<num_channels;i++) {
	previous_index_in_bq[i] = 0;
      }
      for (i=0;i<num_channels;i++) {
	istep_interval[i] = 1;
      }
    }
    /*
      Save the loop state scalars.
    */
    loop_state->num_channels = num_channels;
    loop_state->nx  = nx;
    loop_state->ny  = ny;
    loop_state->nz  = nz;
    loop_state->u_a = u_a;
    loop_state->u_b = u_b;
    loop_state->v_a = v_a;
    loop_state->v_b = v_b;
    loop_state->mu_a = mu_a;
    loop_state->mu_b = mu_b;
    loop_state->mv_a = mv_a;
    loop_state->mv_b = mv_b;
    loop_state->three_mv_a = 3.0 * mv_a;
    loop_state->three_mv_b = 3.0 * mv_b;
    loop_state->recip_mu_a = 1.0/mu_a;
    loop_state->recip_mu_b = 1.0/mu_b;
    loop_state->recip_mv_a = 1.0/mv_a;
    loop_state->recip_mv_b = 1.0/mv_b;
    loop_state->xmin_c = xmin_c;
    loop_state->ymin_c = ymin_c;
    loop_state->xmax_c = xmax_c;
    loop_state->ymax_c = ymax_c;
    loop_state->xc_c   = xc_c;
    loop_state->yc_c   = yc_c;
    loop_state->xd_c   = xd_c;
    loop_state->yd_c   = yd_c;
    loop_state->x_top = x_top;
    loop_state->y_top = y_top;
    loop_state->z_top = z_top;
    loop_state->zmax  = zmax;
    loop_state->lfp   = state->lfp;
    loop_state->max_index = max_index;
    loop_state->min_fract = min_fract;
    loop_state->dvb       = dvb;
    loop_state->dub       = dub;
    loop_state->z_dom_xsum   = z_dom_xsum;
    loop_state->z_dom_ysum   = z_dom_ysum;
    loop_state->z_dom_zsum   = z_dom_zsum;
    loop_state->y_dom_xsum  = y_dom_xsum;
    loop_state->y_dom_ysum  = y_dom_ysum;
    loop_state->y_dom_zsum  = y_dom_zsum;
    loop_state->x_dom_xsum    = x_dom_xsum;
    loop_state->x_dom_ysum    = x_dom_ysum;
    loop_state->x_dom_zsum    = x_dom_zsum;
    loop_state->y_dom_vb_max = y_dom_vb_max;
    loop_state->y_dom_vb_min = y_dom_vb_min;
    loop_state->y_dom_ub_max = y_dom_ub_max;
    loop_state->y_dom_ub_min = y_dom_ub_min;
    loop_state->x_dom_vb_max   = x_dom_vb_max  ;
    loop_state->x_dom_vb_min   = x_dom_vb_min  ;
    loop_state->x_dom_ub_max   = x_dom_ub_max  ;
    loop_state->x_dom_ub_min   = x_dom_ub_min  ;
    loop_state->va_delta     = va_delta;
    loop_state->ua_delta     = ua_delta;
    loop_state->unique_prime_factors = unique_prime_factors;
    loop_state->upf_index    = upf_index;
  }
  return (success);
}      
