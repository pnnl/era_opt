#include "system_includes.h"
#include "state_struct.h"
#include "alloc1.h"
int alloc1(struct state_struct *state) {
  /*
    Allocate the non scalar fields of the state_struct that
    have not yet been allocated:
    Called by: era_opt, gen_combos
    Calls:     calloc, sizeof fprintf,fflush
    
    Uses the following fields of state:
      source_rows,
      source_columns,
      channels_rows,
      channels_columns,
      num_samples,
      num_channels,
      num_huebins,
      num_radiator_temps,
      num_cct_rows,
      subset_size,
      
    Sets the following fields of state:
      channel_matrix,
      channel_matrix_t,
      channel_subset,
      select_counts,      
      ls_counts,
      rs_counts,
      select_max,
      ls_max,
      rs_max,
      steps,
      source_matrix,   
      source_matrix_t, 
      sbar_t_channels, 
      wavelength,  

      scaled_steps,    if (select_only = 0)
      s026_matrix, if (select_only = 0)
      s026,        if (select_only = 0)
      huebins      if (select_only = 0)
      t_a_b_c_s,   if (select_only = 0)
      sbar_10_t_channels if (select_only = 0)
      rx_sbar_10         if (select_only = 0)
      rx_sbar_10_t_channels  if (select_only = 0)
      samples,     if (select_only = 0)
      rx           if (select_only = 0)
      mcat02,      if (select_only = 0)
      mcat02_lu,   if (select_only = 0)
      mhpe,        if (select_only = 0)
      mcat02_ipiv, if (select_only = 0)
      cct_matrix,  if (select_only = 0)
      tt,          if (select_only = 0)
      table_planck, if (select_only = 0)
      mu_wavlength, if (select_only = 0)
      tplanck_num,  if (select_only = 0)
      c2_by_wavelength, if (select_only = 0)
      tplanck,          if (select_only = 0)
      pt_ref_all,       if (select_only = 0)
      pt_test_all,      if (select_only = 0)
      ab_ref_test,      if (select_only = 0)
      ref_hue_bin_angle, if (select_only = 0)
      mean_ab_ref_test,  if (select_only = 0)
      h_bin_ref,         if (select_only = 0)
      xy_norm,           if (select_only = 0)
      gamut_spokes,      if (select_only = 0)
      gamut_edges        if (select_only = 0)
      chroma_shift_rr,   if (select_only = 0)
      mean_ab_diffs,     if (select_only = 0)
      mean_ab_ref_len,   if (select_only = 0)
      ybar_t_channels,   if (select_only = 0)
      ybar_10_t_channels, if (select_only = 0)
      e_t_channels,       if (select_only = 0)
      e_t,                if (select_only = 0)
      sz,                 if (select_only = 0)
      sz_t_channels,      if (select_only = 0)
      bin_count,          if (select_only = 0)
      deltae_samples,     if (select_only = 0)
      deltae_bins,        if (select_only = 0)
      rf_h,               if (select_only = 0)
      href_sum,           if (select_only = 0)
      x_test_norm,        if (select_only = 0)
      y_test_norm,        if (select_only = 0)
      x_ref_norm,         if (select_only = 0)
      y_ref_norm,         if (select_only = 0)
      ref_ucs_table,      if (select_only = 0)
  */
  double *channel_matrix;
  double *channel_matrix_t;
  double *channel_subset;
  double *source_matrix;
  double *source_matrix_t;
  double *cct_matrix;
  double *s026_matrix;
  double *s026;
  double *sz;
  double *sz_t_channels;
  /*
  double *quadrangles;
  */
  double *huebins;
  /*
  double *q_a_b_c;
  double *q_s;
  */
  double *t_a_b_c_s;
  double *sbar_t_channels;
  double *sbar_10_t_channels;
  double *rx_sbar_10;
  double *rx_sbar_10_t_channels;
  double *mcat02;
  double *mcat02_lu;
  double *mcat02_i;
  double *mhpe;
  double *mhpe_mcat02_i;
  double *tt;
  double *wavelength;
  double *table_planck;
  double *planck_rad_work;
  double *mu_wavelength;
  double *tplanck_num;
  double *c2_by_wavelength;
  double *tplanck;
  double *stest;
  double *sref;
  double *steps;
  double *psteps;
  double *dbg_steps;
  double *tm_30_results;
  double *scaled_steps;
  double *pt_ref_all;
  double *pt_test_all;
  double *ab_ref_test;
  double *mean_ab_ref_test;
  double *gamut_spokes;
  double *gamut_edges;
  double *chroma_shift_rr;  /* 2*num_huebins */
  double *mean_ab_diffs;    /* 2*num_huebins */
  double *mean_ab_ref_len;  /* num_huebins */
  double *h_bin_ref;
  double *xy_norm;
  double *ybar_t_channels;
  double *ybar_10_t_channels;
  double *e_t_channels;
  double *e_t;
  double *deltae_samples;
  double *ref_hue_bin_angle; /* num_samples */
  double *deltae_bins;
  double *rf_h;
  double *href_sum;
  double *x_test_norm;
  double *y_test_norm;
  double *x_ref_norm;
  double *y_ref_norm;
  /*
  double *sref_table;
  double *pw_ref_table;
  double *pt_ref_table;
  */
  double *ref_ucs_table;
  double *interp_ucs_row;
  /*
  double *check_ansi_work;
  */
  int64_t *select_counts;
  int64_t *ls_counts;
  int64_t *rs_counts;
  int64_t *select_max;
  int64_t *ls_max;
  int64_t *rs_max;
  
  int    *mcat02_ipiv;
  int    *bin_count;
  int64_t one_l;
  int64_t value_size;
  int64_t ask_for;
  int64_t usage;
  int64_t source_matrix_size;
  int64_t channels_matrix_size;
  int64_t channels_subset_size;
  
  int64_t source_rows;
  int64_t source_columns;

  int64_t channels_rows;
  int64_t channels_columns;

  int64_t s026_columns;

  int64_t num_channels;
  /*
  int64_t num_quadrangles;
  */

  int64_t num_huebins;
  int64_t num_samples;

  int64_t coord_size;
  int64_t eq_size;

  int64_t num_radiator_temps;
  int64_t num_cct_rows;

  int64_t subset_size;

  int success;
  int sz_index;
  

  int rows_ucs_table;
  int t_t_inc;
  int sref_t_t_min;
  int sref_t_t_max;

  int num_labels;
  int select_only;

  
  FILE *lfp;
  FILE *efp;

  success = 1;
  one_l   = (int64_t)1;
  /*
    Unapck state arguments.
  */
  source_rows        = (int64_t)state->source_rows;
  source_columns     = (int64_t)state->source_columns;

  channels_rows      = (int64_t)state->channels_rows;
  channels_columns   = (int64_t)state->channels_columns;

  s026_columns       = (int64_t)state->s026_columns;

  num_channels       = (int64_t)state->num_channels;
  num_samples        = (int64_t)state->num_samples;
  /*
  num_quadrangles    = (int64_t)state->num_quadrangles;
  */
  num_radiator_temps = (int64_t)state->num_radiator_temps;
  num_cct_rows       = (int64_t)state->num_cct_rows;
  
  num_huebins        = (int64_t)state->num_huebins;
  subset_size        = (int64_t)state->subset_size;
  select_only        = state->select_only;
  sz_index           = state->sz_index;
  usage              = state->usage;
  lfp                = state->lfp;
  t_t_inc            = state->t_t_inc;
  sref_t_t_min       = state->sref_t_t_min;
  sref_t_t_max       = state->sref_t_t_max;
  rows_ucs_table    = ((sref_t_t_max - sref_t_t_min)/t_t_inc)+1;
  if (rows_ucs_table < 1) {
    rows_ucs_table = 1;
  }
  state->rows_ucs_table = rows_ucs_table;

  value_size = (int64_t)sizeof(double);
  if (success) {
    channels_matrix_size = channels_rows * channels_columns * value_size;
    ask_for = channels_matrix_size;
    usage += ask_for;
    channel_matrix = (double *)calloc(one_l,ask_for);
    if (channel_matrix == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
		"channel_matrix\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    usage += ask_for;
    channel_matrix_t = (double *)calloc(one_l,ask_for);
    if (channel_matrix_t == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
		"channel_matrix_t\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    channels_subset_size = channels_rows * subset_size * value_size;
    ask_for = channels_subset_size;
    usage += ask_for;
    channel_subset  = (double *)calloc(one_l,ask_for);
    if (channel_subset == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
		"channel_subset\n",ask_for);
	fflush(lfp);
      }
    }
  }    
  if (success) {
    ask_for = ((int64_t)384) * sizeof(int64_t);
    usage += ask_for;
    select_counts = (int64_t*)calloc(one_l,ask_for);
    if (select_counts == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
		"select_counts and select_max\n",ask_for);
	fflush(lfp);
      }
    } else {
      ls_counts  = select_counts +64; /* Caution address arithmetic */
      rs_counts  = ls_counts + 64;    /* Caution address arithmetic */
      select_max = rs_counts + 64;    /* Caution address arithmetic */
      ls_max     = select_max + 64;   /* Caution address arithmetic */
      rs_max     = ls_max + 64;       /* Caution address arithmetic */
    }
  }
  if (success) {
    ask_for = (num_channels + 2 + num_huebins) * value_size;
    usage += ask_for;
    steps = (double*)calloc(one_l,ask_for);
    if (steps == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
		"steps\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    source_matrix_size = source_rows * source_columns * value_size;
    ask_for = source_matrix_size;
    usage   += ask_for;
    source_matrix = (double *)calloc(one_l,ask_for);
    if (source_matrix == NULL) {
      success = 0;
      if (lfp) {
  	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
  		"source_matrix\n",ask_for);
  	fflush(lfp);
      }
    }
  } 
  if (success) {
    usage += ask_for;
    source_matrix_t = (double *)calloc(one_l,ask_for);
    if (source_matrix_t == NULL) {
      success = 0;
      if (lfp) {
  	fprintf(lfp,"thous_opt_2_alloc1: unable to allocate %ld bytes for "
  		"source_matrix_t\n",ask_for);
  	fflush(lfp);
      }
    }
  } 
  if (success) {
    ask_for = 3 * num_channels * value_size;
    usage += ask_for;
    sbar_t_channels = (double *)calloc(one_l,ask_for);
    if (sbar_t_channels == NULL) {
      success = 0;
      if (lfp) {
  	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
  		"sbar_t_channels\n",ask_for);
  	fflush(lfp);
      }
    }
  } 
  if (success) {
    ask_for = source_rows * value_size;
    usage += ask_for;
    wavelength = (double*)calloc(one_l,ask_for);
    if (wavelength == NULL) {
      success = 0;
      if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"tt\n",ask_for);
    	fflush(lfp);
      }
    }
  } 
  /*
    These remaining fields are only used if select_only is 0.
    These are used in the ies_tm_30_18 routine which not called if
    select_only is 1.
  */
  if (select_only == 0) {

    if (success) {
      ask_for = source_rows * s026_columns * value_size;
      usage += ask_for;
      s026_matrix = (double *)calloc(one_l,ask_for);
      if (s026_matrix == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for s026_matrix\n",
    		ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = source_rows * s026_columns * value_size;
      usage += ask_for;
      s026 = (double *)calloc(one_l,ask_for);
      if (s026 == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for s026\n",
    		ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_channels * value_size;
      usage += ask_for;
      sz_t_channels = (double*)calloc(one_l,ask_for);
      if(sz_t_channels == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,
    		"alloc1: unable to allocate %ld bytes for sz_t_channels\n",
    		ask_for);
    	fflush(lfp);
        }
      }
    } 
    /*
    if (success) {
      coord_size = (8 * num_quadrangles);
      eq_size    = (12 * num_quadrangles);
      side_size  = (4  * num_quadrangles);
      ask_for = (coord_size + eq_size + side_size)*value_size;
      usage += ask_for;
      quadrangles = (double*)calloc(one_l,ask_for);
      if (quadrangles == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"quadrangles data\n",ask_for);
    	fflush(lfp);
        }
      } else {
        q_a_b_c = quadrangles + coord_size; // Caution, address arithmetic 
        q_s     = q_a_b_c + eq_size;        // Caution, address arithmetic 
      }
    }
    */
    if (success) {
      coord_size = 6 * num_huebins;
      eq_size    = 12 * num_huebins;
      ask_for = (coord_size + eq_size) * value_size;
      usage += ask_for;
      huebins = (double*)calloc(one_l,ask_for);
      if (huebins == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"huebins data\n",ask_for);
    	fflush(lfp);
        }
      } else {
        t_a_b_c_s = huebins + coord_size; /* Caution, address arithmetic */
      }
    } 
    if (success) {
      ask_for = 3 * num_channels * value_size;
      usage += ask_for;
      sbar_10_t_channels = (double *)calloc(one_l,ask_for);
      if (sbar_10_t_channels == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"sbar_10_t_channels\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = 3 * num_samples * source_rows * value_size;
      usage += ask_for;
      rx_sbar_10 = (double *)calloc(one_l,ask_for);
      if (rx_sbar_10 == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"rx_sbar_10\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for  = 3*num_samples * num_channels * value_size;
      usage += ask_for;
      rx_sbar_10_t_channels = (double *)calloc(one_l,ask_for);
      if (rx_sbar_10_t_channels == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"rx_sbar_10_t_channels\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      /*
        5 3x3 matrices and 1 3x1 vector although with the new method 
        we no longer need mcat02_ipiv.
      */
      ask_for = 48 * value_size;
      usage  +=  ask_for;
      mcat02 = (double *)calloc(one_l,ask_for);
      if (mcat02 == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"mcat02\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    /*
    if (success) {
      ask_for = 4 * num_quadrangles * value_size;
      usage += ask_for;
      check_ansi_work = (double *)calloc(one_l,ask_for);
      if (check_ansi_work == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"thous_opt_2_alloc1: unable to allocate %ld bytes for "
    		"check_ansi_work\n",ask_for);
    	fflush(lfp);
        }
      }
    }
    */
    if (success) {
      mcat02_lu = (double*)&mcat02[9]; 
      mcat02_i  = (double*)&mcat02_lu[9];
      mhpe      = (double*)&mcat02_i[9];
      mhpe_mcat02_i = (double*)&mhpe[9];
      mcat02_ipiv = (int*)&mhpe_mcat02_i[9];
      
      ask_for = num_radiator_temps * value_size;
      usage   += ask_for;
      tt      = (double*)calloc(one_l,ask_for);
      if (tt == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"tt\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_samples * value_size;
      usage += ask_for;
      ref_hue_bin_angle = (double *)calloc(one_l,ask_for);
      if (ref_hue_bin_angle == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"ref_hue_bin_angle\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      /*
        cct_matrix will be a num_cct_rows x 3 matrix stored in 
        row majort order.
      */
      ask_for = num_cct_rows * 3 * value_size;
      usage += ask_for;
      cct_matrix = (double *)calloc(one_l,ask_for);
      if (cct_matrix == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"cct_matrix\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = 4 * num_radiator_temps * value_size;
      usage   += ask_for;
      table_planck = (double*)calloc(one_l,ask_for);
      if (table_planck == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"table_planck\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      /*
        planck_rad_work, we need mu_wavelength = wavelength * .0000000001
                         tplanck_num = c1 ./ (mu_wavelength .^ 5)
    		       c2_by_wavelength = c2 ./mu_wavelength
    		       tplanck
      */
      ask_for = source_rows * 4 * value_size;
      usage += ask_for;
      planck_rad_work =(double*)calloc(one_l,ask_for);
      if (planck_rad_work == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"planck_rad_work\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    /*
      Allocate space for additional vectors use in ies_tm_30_18
    */
    if (success) {
      mu_wavelength    = planck_rad_work;
      /* Caution the next three statments do address arithmetic; */
      tplanck_num      = mu_wavelength + source_rows; 
      c2_by_wavelength = tplanck_num + source_rows;
      tplanck          = c2_by_wavelength + source_rows;
      ask_for = source_rows * value_size;
      usage += ask_for;
      stest = (double*)calloc(one_l,ask_for);
      if (stest == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"stest\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = source_rows * value_size;
      usage += ask_for;
      sref = (double*)calloc(one_l,ask_for);
      if (sref == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"sref\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = (2 + num_huebins) * value_size;
      usage += ask_for;
      tm_30_results = (double*)calloc(one_l,ask_for);
      if (tm_30_results == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"tm_30_results\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_channels * value_size;
      usage += ask_for;
      scaled_steps = (double*)calloc(one_l,ask_for);
      if (scaled_steps == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"scaled_steps\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = 3*num_samples * value_size;
      usage += ask_for;
      pt_ref_all = (double*)calloc(one_l,ask_for);
      if (pt_ref_all == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"pt_ref_all\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = 3*num_samples * value_size;
      usage += ask_for;
      pt_test_all = (double*)calloc(one_l,ask_for);
      if (pt_test_all == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"pt_test_all\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = 4*num_samples * value_size;
      usage += ask_for;
      ab_ref_test = (double*)calloc(one_l,ask_for);
      if (ab_ref_test == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"ab_ref_test\n",ask_for);
    	fflush(lfp);
        }
      }
    }       
    if (success) {
      ask_for = 4*num_huebins * value_size;
      usage += ask_for;
      mean_ab_ref_test = (double*)calloc(one_l,ask_for);
      if (mean_ab_ref_test == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"mean_ab_ref_test\n",ask_for);
    	fflush(lfp);
        }
      }
    }       
    if (success) {
      ask_for = 2*num_huebins * value_size;
      usage += ask_for;
      gamut_spokes = (double*)calloc(one_l,ask_for);
      if (gamut_spokes == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"gamut_spokes\n",ask_for);
    	fflush(lfp);
        }
      }
    }       
    if (success) {
      ask_for = 2*num_huebins * value_size;
      usage += ask_for;
      gamut_edges = (double*)calloc(one_l,ask_for);
      if (gamut_edges == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"gamut_edges\n",ask_for);
    	fflush(lfp);
        }
      }
    }       
    if (success) {
      ask_for = 2*num_huebins * value_size;
      usage += ask_for;
      chroma_shift_rr = (double*)calloc(one_l,ask_for);
      if (chroma_shift_rr == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"chroma_shift_rr\n",ask_for);
    	fflush(lfp);
        }
      }
    }       
    if (success) {
      ask_for = 2*num_huebins * value_size;
      usage += ask_for;
      mean_ab_diffs= (double*)calloc(one_l,ask_for);
      if (mean_ab_diffs == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"mean_ab_diffs\n",ask_for);
    	fflush(lfp);
        }
      }
    }       
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      mean_ab_ref_len= (double*)calloc(one_l,ask_for);
      if (mean_ab_ref_len == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"mean_ab_ref_Len\n",ask_for);
    	fflush(lfp);
        }
      }
    }       
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      h_bin_ref = (double*)calloc(one_l,ask_for);
      if (h_bin_ref == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"h_bin_ref\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = 4*num_huebins*value_size;
      usage += ask_for;
      xy_norm = (double*)calloc(one_l,ask_for);
      if (xy_norm == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"xy_norm\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_channels *value_size;
      usage += ask_for;
      ybar_t_channels = (double*)calloc(one_l,ask_for);
      if (ybar_t_channels == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"ybar_t_channels\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_channels *value_size;
      usage += ask_for;
      ybar_10_t_channels = (double*)calloc(one_l,ask_for);
      if (ybar_10_t_channels == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"ybar_10_t_channels\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_channels *value_size;
      usage += ask_for;
      e_t_channels = (double*)calloc(one_l,ask_for);
      if (e_t_channels == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"e_t_channels\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = channels_rows*value_size;
      usage += ask_for;
      e_t = (double*)calloc(one_l,ask_for);
      if (e_t == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"e_t\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_samples * value_size;
      usage += ask_for;
      deltae_samples = (double*)calloc(one_l,ask_for);
      if (deltae_samples == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"deltae_samples\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      deltae_bins = (double*)calloc(one_l,ask_for);
      if (deltae_bins == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"deltae_bins\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      rf_h = (double*)calloc(one_l,ask_for);
      if (rf_h == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"rf_h\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      href_sum = (double*)calloc(one_l,ask_for);
      if (rf_h == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"href_sum\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      x_ref_norm = (double*)calloc(one_l,ask_for);
      if (rf_h == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"x_ref_norm\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      y_ref_norm = (double*)calloc(one_l,ask_for);
      if (rf_h == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"y_ref_norm\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      x_test_norm = (double*)calloc(one_l,ask_for);
      if (rf_h == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"x_test_norm\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = num_huebins * value_size;
      usage += ask_for;
      y_test_norm = (double*)calloc(one_l,ask_for);
      if (rf_h == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"y_test_norm\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    /*
    if (success) {
      ask_for = rows_ucs_table * source_rows * value_size;
      usage += ask_for;
      sref_table = (double*)calloc(one_l,ask_for);
      if (sref_table == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"sref_table\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = rows_ucs_table * 3 * num_samples * value_size;
      usage += ask_for;
      pt_ref_table = (double *)calloc(one_l,ask_for);
      if (pt_ref_table == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"pt_ref_table\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
     if (success) {
       ask_for = rows_ucs_table * 3 * value_size;
       usage += ask_for;
       pw_ref_table = (double *)calloc(one_l,ask_for);
       if (pw_ref_table == NULL) {
         success = 0;
         if (lfp) {
     	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
     		"pw_ref_table\n",ask_for);
     	fflush(lfp);
         }
       }
     }
    */ 
    if (success) {
      ask_for = rows_ucs_table * 5 * num_samples * value_size;
      usage += ask_for;
      ref_ucs_table = (double *)calloc(one_l,ask_for);
      if (ref_ucs_table == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"ref_ucs_table\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    if (success) {
      ask_for = 5 * num_samples * value_size;
      usage += ask_for;
      interp_ucs_row = (double *)calloc(one_l,ask_for);
      if (interp_ucs_row == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"interp_ucs_row\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
      
    if (success) {
      ask_for = num_huebins * sizeof(int);
      usage += ask_for;
      bin_count = (int*)calloc(one_l,ask_for);
      if (bin_count == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"bin_count\n",ask_for);
    	fflush(lfp);
        }
      }
    } 
    /*
      isteps_per_channel is allocated in alloc0.
    if (success) {
      ask_for = num_channels * sizeof(int);
      usage += ask_for;
      isteps_per_channel = (int*)calloc(one_l,ask_for);
      if (isteps_per_channel == NULL) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"alloc1: unable to allocate %ld bytes for "
    		"isteps_per_channel\n",ask_for);
    	fflush(lfp);
        }
      }
    }
    */
  }
  if (success) {
    state->channel_matrix    	 = channel_matrix;
    state->channel_matrix_t  	 = channel_matrix_t;
    state->channel_subset        = channel_subset;
    state->steps              	 = steps;
    state->select_counts         = select_counts;
    state->ls_counts             = ls_counts;
    state->rs_counts             = rs_counts;
    state->select_max            = select_max;
    state->ls_max                = ls_max;
    state->rs_max                = rs_max;
    state->source_matrix      	 = source_matrix;
    state->source_matrix_t    	 = source_matrix_t;
    state->sbar_t_channels    	 = sbar_t_channels;
    state->wavelength         	 = wavelength;

    if (select_only == 0) {
      state->samples            	 = source_matrix_t;
      state->rx                 	 = source_matrix_t;
      /*
      state->quadrangles        	 = quadrangles;
      */
      state->huebins            	 = huebins;
      state->t_a_b_c_s           	 = t_a_b_c_s;
      state->sbar_10_t_channels 	 = sbar_10_t_channels;
      state->rx_sbar_10         	 = rx_sbar_10;
      state->rx_sbar_10_t_channels = rx_sbar_10_t_channels;
      /*
      state->check_ansi_work    	 = check_ansi_work;
      */
      state->mcat02             	 = mcat02;
      state->mcat02_lu           	 = mcat02_lu;
      state->mcat02_i                    = mcat02_i;
      state->mhpe               	 = mhpe;
      state->mhpe_mcat02_i               = mhpe_mcat02_i;
      state->mcat02_ipiv        	 = mcat02_ipiv;
      state->cct_matrix                  = cct_matrix;
      state->tt                 	 = tt;
      state->table_planck       	 = table_planck;
      state->mu_wavelength      	 = mu_wavelength;
      state->tplanck_num        	 = tplanck_num;
      state->c2_by_wavelength   	 = c2_by_wavelength;
      state->tplanck            	 = tplanck;
      state->stest              	 = stest;
      state->sref               	 = sref;
      state->ybar_t_channels       = ybar_t_channels;
      state->ybar_10_t_channels    = ybar_10_t_channels;
      state->e_t_channels          = e_t_channels;
      state->e_t                   = e_t;
      state->s026_matrix           = s026_matrix;
      state->s026                  = s026;
      state->sz_t_channels         = sz_t_channels;                   
      /*
#define DBG 1
      */
#ifdef DBG
      if (lfp) {
        fprintf(lfp,"alloc1: state->sref = %p\n",state->sref);
        fflush(lfp);
      }
#endif 
      state->tm_30_results         = tm_30_results;
      state->scaled_steps       	 = scaled_steps;
      state->pt_ref_all            = pt_ref_all;
      state->pt_test_all           = pt_test_all;
      state->ab_ref_test           = ab_ref_test;
      state->ref_hue_bin_angle     = ref_hue_bin_angle;
      state->mean_ab_ref_test      = mean_ab_ref_test;
      state->gamut_spokes          = gamut_spokes;
      state->gamut_edges           = gamut_edges;
      state->chroma_shift_rr       = chroma_shift_rr;
      state->mean_ab_diffs         = mean_ab_diffs;
      state->mean_ab_ref_len       = mean_ab_ref_len;
      state->h_bin_ref             = h_bin_ref;
      state->xy_norm               = xy_norm;
      state->deltae_samples        = deltae_samples;
      state->deltae_bins           = deltae_bins;
      state->rf_h                  = rf_h;
      state->href_sum              = href_sum;
      state->x_ref_norm            = x_ref_norm;
      state->y_ref_norm            = y_ref_norm;
      state->x_test_norm           = x_test_norm;
      state->y_test_norm           = y_test_norm;
      state->bin_count             = bin_count;
      /*
      state->sref_table            = sref_table;
      state->pw_ref_table          = pw_ref_table;
      state->pt_ref_table          = pt_ref_table;
      */
      state->ref_ucs_table         = ref_ucs_table;
      state->interp_ucs_row        = interp_ucs_row;
    }
    state->usage              	 = usage;
  }
  return(success);
}
