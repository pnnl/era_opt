#ifndef _STATE_STRUCT_H_
#define _STATE_STRUCT_H_ 2
#include "mpi.h"
#include "loop_state_struct.h"
/*
  This is now included in groups_struct.
#include "ellipse_struct.h"
*/
#include "groups_struct.h"
#include "comm_map_struct.h"
struct state_struct {
  struct loop_state_struct *loop_state;
  struct groups_struct      *groups_state;
  struct comm_map_struct    *reduction_map;
  char *parameter_file;      /* 128 allocated in alloc0 */
  char *log_file;            /* 128 allocated in alloc0 */
  char *channels_file;       /* 128 allocated in alloc0 */
  char *subset_file;         /* 128 allocated in alloc0 */
  char *z_dom_file;           /* 128 allocated in alloc0 */
  char *y_dom_file;          /* 128 allocated in alloc0 */
  char *x_dom_file;          /* 128 allocated in alloc0 */
  char *combos_out_file;     /* 128 allocated in alloc0 */
  char *output_file;         /* 128 allocated in alloc0 */
  char *output_dir;          /* 128 allocated in alloc0 */
  char *combos_in_file;      /* 128 allocated in alloc0 */
  char *source_file;         /* 128 allocated in alloc0 */
  char *cct_file;            /* 128 allocated in alloc0 */
  /*
  char *quadrangles_file;    
  */
  char *huebins_file;        /* 128 allocated in alloc0 */
  char *s026_file;           /* 128 allocated in alloc0 */
  char *summary_file;        /* 128 allocated in alloc0 */
  char *input_dir;           /* 128 allocated in alloc0 */
  char *combo_dir;           /* 128 allocated in alloc0 */
  char *mcat02inv_file;      /* 128 allocated in alloc0 */

  char *dir_plus_fn;         /* 256 allocated in alloc0 */

  char *io_buffer;           /* 16384 allocated in alloc0 */
  double *source_matrix;   /* source_rows x source_columns stored in row major order */
  double *channel_matrix;  /* channels_rows x channels_columns (row major) */
  double *channel_subset;  /* Work space for column subset of the channels matrix, 
			      channels_rows * subset_size */
  double *full_channels;   /* pointer to full channels matrix so that channel_matrix may be used to
			      point to a channel_subbset */
  double *s026_matrix;      /* source_rows x 8 in row major order */
  double *source_matrix_t; /* source_rows x source_columns (column major order) */
  double *channel_matrix_t; /* channels_columns x channels_rows */
  double *channels;      /* channels_rows x channels columns in column major order */
  double *all_channels;  /* channels_rows x channels columns in column major order */
                         /* all_channels is a pointer to the full set of channels */
                         /* when a subset is being used */
  double *s026;           /* source_rows x 8 in column major order */  
			    
  /*double *quadrangles; */
  double *huebins;   /* 6*num_huebins */
  /*
  double *q_a_b_c; // 12 x num_quadrangles 
  double *q_a;
  double *q_b;
  double *q_c;
  double *q_s;   // 4 x num_quadrangles 
  */
  double *cct_matrix;   /* num_cct_rows x 3 */
  double *tt; /* radiator temperatures, length is 951 1000:20:10000 
		 Allocated in alloc1, initialized in prep,
		 used in planck_rad, */
  double *wavelength; /* wavelengths vector, length 401  380:1:780 
		 Allocated in alloc1, initialized in prep,
		 used in planck_rad, */
  double *table_planck; /* tt,ui,vi,di_sq table.
		 Allocated in alloc1, initialized in prep,
		 used in planck_rad, */
  double *t_a_b_c_s; /* 12 x num_huebins */
  double *sbar_t_channels; /* 3 x num_channels  = (xbar,ybar,zbar)' * channels 
			    stored in row major order */

  double *ybar_t_channels; /* num_channels, extracted from sbar_t_channels,
			      used in computing LER value. */

  double *ybar_10_t_channels; /* num_channels, extracted from sbar_10_t_channels,
			      used in computing stest_scale value. */

  double *e_t_channels;    /* num_channels, sum of the columns of the channels
			      matrix, used in computing LER value. */

  double *e_t;             /* channel_rows long vector of ones. */
  double *sz;              /* Melanopic vector from s026 (fourth column)
			      Length is source_rows */
  double *sz_t_channels;   /* Length is num channels*/

  double *sbar_10_t_channels; /* 3 x num_channels  = 
			      (xbar_10,ybar_10,zbar_10)' * channels 
			    stored in row major order */
  double *rx_sbar_10_t_channels; /* 3*num_samples x num_channels */
  double *sbar;    /* xbar, ybar, zbar rows of source matrix^T */
  double *sbar_10; /* xbar_10, ybar_10, zbar_10 rows of source matrix^T */
  double *rx_sbar_10; /* source_rows x 3 * num_samples, stored by columns */
  double *xbar;    /* source_rows */
  double *ybar;    /* source_rows */
  double *zbar;    /* source_rows */
  double *xbar_10; /* source_rows */
  double *ybar_10; /* source_rows */
  double *zbar_10; /* source_rows */
  double *s_0;     /* source_rows */
  double *s_1;     /* source_rows */
  double *s_2;     /* source_rows */
  double *samples;  /* source_file_rows x num_samples stored in column major order */
  double *rx;       /* synonum for samples */
  double *mcat02;  /* mcat02 3 x 3 matrix stored in column major order. */
  double *mcat02_lu; /* lu factorization of mcat02, length is 9 */
  double *mcat02_i;  /* Inverse of mcat02 used in computing mhpe mcat02^-1 product. */
  double *mhpe;     /* mhpe 3 x 3 matrix stored in column major order. */
  double *mhpe_mcat02_i; /* mhpe * mcat02_i used in cam02ucs */
  /*
      planck_rad_work, we need mu_wavelength = wavelength * .0000000001
                       tplanck_num = c1 ./ (mu_wavelength .^ 5)
		       c2_by_wavelength = c2 ./mu_wavelength
		       tplanck
  */
  double *mu_wavelength;    /* source_rows */
  double *tplanck_num;      /* source_rows */
  double *c2_by_wavelength; /* source_rows */
  double *tplanck;          /* source_rows */
  double *stest;            /* source_rows */
  double *sref;             /* source_rows */
  double *steps;            /* num_channels */
  double *dbg_steps;        /* num_channels */
  double *tm_30_results;    /* num_huebins+2: rf, rg, r_c_sh(0:num_huebins-1) */
  double *scaled_steps;     /* num_channels */
  double *pt_ref_all;       /* 3*num_samples */
  double *pt_test_all;      /* 3*num_samples */
  double *ab_ref_test;      /* 4*num_samples */
  double *ref_hue_bin_angle; /* num_samples */
  double *mean_ab_ref_test; /* 4*num_huebins */
  double *gamut_ref_test;   /* synonym for mean_ab_ref_test */
  double *gamut_spokes;     /* 2*num_huebins */
  double *gamut_edges;      /* 2*num_huebins */
  double *chroma_shift_rr;  /* 2*num_huebins */
  double *mean_ab_diffs;    /* 2*num_huebins */
  double *mean_ab_ref_len;  /* num_huebins */
  double *h_bin_ref;        /* num_huebins Used to sum the arc-tangents (in 
			       degress) of ab_ref sample points in ies_tm_30 */
  double *deltae_bins;      /* num_huebins: sum of the cam02ucs_deltaes landing in a
			       hue bin, and later their average */
  double *deltae_samples;   /* num_samples: stores the cam02ucs_deltae 
			       for each sample */
  double *rf_h;             /* num_huebins: Local color fidelity computed 
			       from deltae_bins */
  double *xy_norm;          /* 4*num_huebins */
  double *shifted_x;        /* num_huebins */
  double *shifted_y;        /* num_huebins */
  double *href_sum;         /* num_huebins */
  double *x_ref_norm;       /* num_huebins */
  double *y_ref_norm;       /* num_huebins */
  double *x_test_norm;      /* num_huebins */
  double *y_test_norm;      /* num_huebins */
  double *scales;           /* number_scales (refinement levels) */

  double *vrb_radii;    /* Test radii for number of channels in 
			   vrb_first_indes, vrb_last index range */

  double *test_radii_sq;   /* Length is number channels, square of test radius
			      for that number of channels, set from vrb_* 
			      arrays. */
  double *x_dom_data;        /* Pointer to x_dom (red) combos in mmaped x_dom_file */
  double *y_dom_data;      /* Pointer to y_dom (green) combos in mmaped y_dom_file */
  double *z_dom_data;       /* Pointer to z_dom (blue) combos in mmaped z_dom_file */

  
  //double *sref_table;   /* Table of sref vectors  for lookup or interpolation 
  //			     source_rows x rows_ucs_table
  //		          */
  //double *pw_ref_table; /* Table of pw_ref vectors 3 x rows_ucs_table */
  //double *pt_ref_table; /* Table of pt_ref vectors 3*num_samples x rows_ucs_table */
  
  /* 
     Replace the above three with ref_ucs_table to save a table of the for ucs outputs
     from the cam02ucs call with sref and each of the source sample SPDs.
     Size is 5*num_samples * rows_ucs_table
     the 5 ucs_coords per sample are japos,mapos,aaposm,baposm, and hue_bin_angle
  */
  double *ref_ucs_table;
  double *interp_ucs_row; /* Length is 5*num_samples used to interpolate between two
			     adjacent rows of the ref_ucs_table leaving it unchanged. */


  int64_t *iy_dom_bins; /* pointer to y-dominant combo bin delimiters in
			   mmaped y_dom_file */
  int64_t *ix_dom_bins;   /* pointer to x-dominant combo bin delimieters in
			   mmaped x_dom_file */

  int64_t *select_counts; /* length is 4 for accumulated ix_delta_sum, iy_delta_sum and num_in_tca + spare */

  int64_t *ls_counts;
  int64_t *rs_counts;
  
  int64_t *select_max;   /* length is 4 for max lix_delta_sum, iy_delta and num_in_tca + spare */
  int64_t *ls_max;
  int64_t *rs_max;

  int *mcat02_ipiv; /* ipiv vector used in forming lu factoriaction of mcat02
		       length is 3. */
  int *bin_count; /* bin_count used to count ref points that land in huebins.
		     lengh is num_huebins. allocated in alloc1. */
  int *isteps_per_channel; /* max_num_channels (64) allocated in alloc0 */

  int *vrb_first_index; /* lowest number of channels test radius applies to */
  
  int *vrb_last_index;  /* Largest number of channels test radius applies to */

  int *subset_indices;  /* Vector of length max_subset_size with indices in [1..num_channels] of subset
		           of indices of channels matrix columns to use as channels matrix with
	                   first_channel_index being counted as index 1. */		      
  /*
    Constants used to define "soft white light" filter.
  */
  double u_target;
  double v_target;
  double test_rad;
  double rf_thresh;
  double rad_skosh;
  /*
    Planckian radiation constants set in prep
  */
  double c1;
  double c2;
  /*
    Constants used in cam02ucs, set in prep
  */
  double pi;
  double two_pi;
  double half_pi;
  double degrees_per_radian;
  double radians_per_degree;
  double n;
  double yb;
  double yw;
  double fl;
  double fl_over_100;
  double fourth_root_fl;
  double nbb;
  double ncb;
  double zz;
  double c;
  double c_zz;
  double f;
  double sr;
  double nc;
  double chroma_c_scale;
  double degree_adaptation;
  double twentieth;
  double eleventh;
  double ninth;
  double twelve_elevenths;
  double hr_thresh;
  double t_const;
  double sin_2;
  double cos_2;

  double y_dom_combo_bin_width;
  double y_dom_vb_min;
  double y_dom_vb_max;
  double y_dom_ub_min;
  double y_dom_ub_max;

  double x_dom_combo_bin_width;
  double x_dom_vb_min;
  double x_dom_vb_max;
  double x_dom_ub_min;
  double x_dom_ub_max;
  double phi_degrees_cor;
  /*
  double svd_phi_degrees;
  */
  double neutral_rcs_radius;  /* if all rcsh huebin lengths are < neutral_rcs_radius, the LCS group is neutral */
  double neutral_skew_radius; /* if major axis length/minor axis_length < 1 + neutral_skew_radius ellipse 
				 group is neutral */

  double watchme; /* frequency of progress report to log file, default value = .05 report every 5 */

  /*
    Scalars to track min cct t_t, max cct t_t, and max rf
  */
  double cct_t_t_min;
  double cct_t_t_max;
  double rf_max;

  int64_t usage;
  int64_t num_in_tca;
  int64_t iy_delta_sum;
  int64_t ix_delta_sum;
  
  int64_t max_in_tca;
  int64_t max_iy_delta_sum;
  int64_t max_ix_delta_sum;
  int64_t max_number_results; /* per node */
  int64_t total_number_results;

  int64_t izc;
  int64_t iyc;

  int64_t ixc;

  int my_id;
  int thread_id;

  int num_procs;
  int num_threads;

  int source_rows;
  int source_columns;

  int channels_rows;
  int channels_columns;

  int first_channel_index;
  int last_channel_index;

  int num_channels;
  int num_huebins;

  int subset_size;
  int max_subset_size; /* declared length of subset_indices vector above 64 */
  
  int buffer_length;
  int num_samples;

  int wavelength_index;
  int xbar_index;

  int ybar_index;
  int zbar_index;

  int xbar_10_index;
  int ybar_10_index;

  int zbar_10_index;
  int s_0_index;

  int s_1_index;
  int s_2_index;

  int sz_index;
  int s026_columns;

  int num_radiator_temps;
  int max_num_channels; /* 64*/

  int nc_steps_set; /* Number of steps per channel set in the
		       STEPS_PER_CHANNEL line of the parameters file.
		    */
  int num_cct_rows;

  int use_multiscale; /* 1 is yes, 0 is no, default is 0 */
  int number_scales; /* Default 1, no refinement */
                            
  int max_scales;
  int max_index;

  int max_radii_bins;
  int extra_log_info;

  int number_vrb; /* Number variable radii bins */
  int use_vrb;    /* 0 for use fixed test_rad radius, 1 for use 
		     variable radius bins */



  /* File descriptors for the memory mapped combo files. */
  int xfd;
  int yfd;

  int zfd;
  int nz_dom_fields;     /* Number of 8byte fields in a z_dom combo record
			   in the mmapped z_dom_file */

  int iz_dom_tsv_offset; /* Tristimulas value offset (in 8 byte words) in
			   a z_dom combo record in the mmapped z_dom_file */
  int iz_dom_sbpi;       /* Number of stored bits per index in the 
			   istep field of a z_dom combo record 
			   in the mmapped z_dom_file, needed for
			   unpacking istep_index values from the combo record */

  int ny_dom_fields;    /* Number of 8byte fields in a y_dom combo record
			   in the mmapped y_dom_file */
  int ny_dom_bins;      /* Number of bins for y_dom_combos sorted by vb_value. */

  int iy_dom_tsv_offset; /* Tristimulas value offset (in 8 byte words) in
			   a y_dom combo record in the mmapped y_dom_file */
  int iy_dom_sbpi;       /* Number of stored bits per index in the 
			   istep field of a green combo record 
			   in the mmapped y_dom_file, needed for
			   upacking istep_index values from the combo record */

  int nx_dom_fields;       /* Number of 8byte fields in a x_dom combo record
			   in the mmapped x_dom_file */
  int nx_dom_bins;         /* Number of bins for x_dom combos sorted by ub_value. */

  int ix_dom_tsv_offset;   /* Tristimulas value offset (in 8 byte words) in
			   a x_dom combo record in the mmapped x_dom_file */
  int ix_dom_sbpi;         /* Number of stored bits per index in the 
			   istep field of a x_dom combo record 
			   in the mmapped x_dom_file, needed for
			   upacking istep_index values from the combo record */

  int max_num_bins;
  int nz_cpcg_max;       /*
			   Maximum number of NonZero Channels 
                           Per Coordinate Group, -1 is use default
			   value = num_channels for no maximum
			 */
  int full_num_channels; /* last_channel_index - first_channel_index + 1 */
  int num_labels; /* currently 16 */

  int fit_ellipse_status;
  int use_bruteforce;

  int filename_length;
  int full_output; /*     full_output, 1 = output_tm_30_results + group summary
			               0 = only group_summary default will be 0.
		   */

  int select_only; /* 1 for don't call ies_tm_30 routines, only count number combos in tca. 
		      0 for call ies_tm_30 routines for combos in tca. Default is 0 */
  int discount_illuminant; /* 1 to use degree_adaptation = 1; 0 to use 
			      D=F*(1-(1/3.6)*exp((-LA-42)/92));
			      Fixed at 1 now in prep.
			   */

  int new_cam02ucs;  /* 1 to use new cam02ucs, 0 to use old one
			default is 0 */

  int print_efit;    /* 1 to cause fit_ellipse to print out intermediates */
  int hp_out;        /* 0 for f6.2 output format, 1 for f16.12 format    */
                     /* should only be used for small number of channels  
		        to avoid prohibitively large volume of output     */

  int save_steps;    /* 1 to save combo steps to combos_file when select_only is chosen */
                     /* default is 0. */
  int read_combos;   /* 1 to read combos in instead of executing era_opt, for testing
			ies_tm_30_18 routine */

  int use_alt_cct; /* 1 for alternate cct calculation */
  int read_mcat02inv; /* 0 for comput mcat02^-1, 
			 1 (default) for read it in from mcat02inv_file */

  int cam02ucs_1mvp;
  int cam02ucs_2mvp;

  int sref_mode;    /* 0 for compute sref and ucs_row
		       from t_t (original)
		       1 for lookup sref, ucs_row 
		       2 for interpolate sref, ucs_row 
		    */
  int t_t_inc;      /* increment value to use in building tables
		       for sref_mod = 1 or 2 */

  int sref_t_t_min; /* minimum t_t value for sref_mode =1 or 2*/
  int sref_t_t_max; /* maximum t_t value for sref_mode =1 or 2*/

  int rows_ucs_table; /* Number of rows in the ref_ucs_table 
	    	         = (sref_t_t_max - sref_t_t_min)/t_t_inc + 1
		      */
  int dump_combos;  /* when select_only is 1, if dump_combos is also 1,
		       then write the combos_out file(s) with
		       the power steps and tsv values for each combination
		       that lands in the TCA. This will be a csv file,
		       the first line just have number_of_channels,number_of_power_levels
		       Subsequent lines will have number_of_channels comma separated integers
		       representing the channel power levels 
		       (these will be integers in the range [0,number_of_power_levels-1].
		       followed by the the three comma separated tsv coordinates (x,y,z)
		    */


  FILE *pfp; /* parameter file */
  FILE *lfp; /* log file */

  FILE *sfp; /* summary file */
  FILE *cfp; /* Combos input file */

  /*
  FILE *qfp; 
  FILE *tfp; 
  */

  FILE *ofp; /* output file */
  FILE *xfp; /* x-dominant combos output file */

  FILE *yfp; /* y-dominant combos output file */
  FILE *zfp; /* z-dominant combos output file */

  FILE *dfp; /* Combos output file. */
  FILE *efp; /* place holder */
}
;  
#endif
