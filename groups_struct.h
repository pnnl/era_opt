#ifndef _GROUPS_STRUCT_H_
#define _GROUPS_STRUCT_H_ 1
#include "ellipse_struct.h"
struct groups_struct {
  /*
    Struct to hold histograms and group counts. and ellipse structs in alloc2 or alloc3.
  */
  struct ellipse_struct     *ellipse;
  /*
  struct ellipse_struct     *svd_ellipse;
  */
  int64_t usage; /* about 90Mbytes */
  double *ellipse_work;   /* 14*huebins */
  double  *pvf_range_min; /* 4 x 9 entries */
  double  *pvf_range_max; /* 4 x 9 entries */
  double  *pvf_range_left; /* 36 */
  double  *pvf_range_right; /* 36 */
 
  int64_t *class1_histograms;    /* (2^18=262144) x num_class1_groups; */
  int64_t *class1_sh;           /* Class 1 subhistograms pvf group counts group index */
  int64_t *class1_unique_counts; /* num_class1_groups */
  int64_t *class2_histograms;    /* 64 x 2 : counts of number of entries in each unique group */
  int64_t *class2_unique_counts;  /* num_class2_groups */
  int64_t *all_counts;           /* length is 32, overlapps p_counts, v_counts and f_counts, */
  int64_t *p_counts;             /* length is 4 */
  int64_t *v_counts;             /* length is 4 */
  int64_t *f_counts;             /* length is 4 */

  /*
  int *class1_groups;            // list of unique labels for each class 1group 
			         //   (2^20=1048576 x num_class1_groups; 
  int *class2_groups;            // 64 x 2 : list of unique PVF labels within in each group. 
  */
  /* Class 1 groups */
  /*
  int *lcs_17_u_d;
  int *lcs_9_u_c;
  int *lcs_9_u_d;
  int *lcs_5_u_c;
  int *lcs_5_u_o;
  int *lcs_9_b_d;
  int *lcs_5_b_c;
  int *lcs_5_b_d;
  int *ell_9_b_d;
  int *ell_5_b_c;
  int *ell_5_b_d;
  */
  /*
  int *svd_9_b_d;
  int *svd_5_b_c;
  int *svd_5_b_d;
  */

  /* Class 1 histograms */
  int64_t *lcs_17_u_d_h;
  int64_t *lcs_9_u_c_h;
  int64_t *lcs_9_u_d_h;
  int64_t *lcs_5_u_c_h;
  int64_t *lcs_5_u_o_h;
  int64_t *lcs_9_b_d_h;
  int64_t *lcs_5_b_c_h;
  int64_t *lcs_5_b_d_h;
  int64_t *ell_9_b_d_h;
  int64_t *ell_5_b_c_h;
  int64_t *ell_5_b_d_h;
  /* Class 1 subhistograms */
  int64_t *lcs_17_u_d_sh;
  int64_t *lcs_9_u_c_sh;
  int64_t *lcs_9_u_d_sh;
  int64_t *lcs_5_u_c_sh;
  int64_t *lcs_5_u_o_sh;
  int64_t *lcs_9_b_d_sh;
  int64_t *lcs_5_b_c_sh;
  int64_t *lcs_5_b_d_sh;
  int64_t *ell_9_b_d_sh;
  int64_t *ell_5_b_c_sh;
  int64_t *ell_5_b_d_sh;
  /*
  int64_t *svd_9_b_d_h;
  int64_t *svd_5_b_c_h;
  int64_t *svd_5_b_d_h;
  */
  int64_t *left_c1_h;
  int64_t *right_c1_h;
  int64_t *left_c1_sh;
  int64_t *right_c1_sh;
  /* Class 2 groups */
  /*
  int *pvfp;
  int *pvfa;
  */
  /* Class 2 histograms */
  int64_t *pvfp_h;
  int64_t *pvfa_h;
  int64_t *left_c2_h;
  int64_t *right_c2_h;

  /* Class 2 subhistograms */
  int64_t *p1_sh;
  int64_t *p2_sh;
  int64_t *p3_sh;
  int64_t *v1_sh;
  int64_t *v2_sh;
  int64_t *v3_sh;
  int64_t *f1_sh;
  int64_t *f2_sh;
  int64_t *f3_sh;
  int64_t *left_p_sh;
  int64_t *right_p_sh;
  int64_t *left_v_sh;
  int64_t *right_v_sh;
  int64_t *left_f_sh;
  int64_t *right_f_sh;

  /* maps of huebin index to group number */
  int *map_9_u_c; /* [17] = {0,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,1} */
  int *map_9_u_d; /* [17] = {0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8} */
  int *map_5_u_c; /* [17] = {0,1,1,2,2,2,2,3,3,3,3,4,4,4,4,1,1} */
  int *map_5_u_o; /* [17] = {0,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,1} */
  int *map_9_b_d; /* [17] = {0,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8} */
  int *map_5_b_c; /* [17] = {0,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4,1} */
  int *map_5_b_d; /* [17] = {0,1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4} */
  int *labels;    /* [16] group labels for: nint_rf, nint_rg,
		     17-U-D, 9-U-C, 9-U-D, 5-U-C, 5-U-O, 
		     9-B-D 5-B-C 5-B-D, P, V, and F */

  /* Note that 17-U-D label is just max huebin index. */

  /* Count accumulation buffers to accumulate the following 5 quantities in addition to the 12 pvf counts
     lcs_neutral_count, ellipse_neutral_count, pvf_neutral_count, number_results
     will allocate each as array of length 32 integers for alignment reasons.
  */
  /*
    svd_neutral_count, not used any more
  */
  int64_t *left_counts;
  int64_t *right_counts;

  int64_t lcs_neutral_count;
  int64_t ellipse_neutral_count;

  int64_t pvf_neutral_count;
  /*
  int64_t svd_neutral_count;
  */
  int64_t number_results;
  int64_t num_in_tca;

  int64_t ix_delta_sum;
  int64_t iy_delta_sum;



  /* scalars */
  int num_class1_groups;
  int length_class1_group; 

  int num_class2_groups;
  int length_class2_group;

  int nint_r_f_bits;   /* number of bits in nint_r_f in [60:100] so 7 bits (0:127)
			  we will subtract 60 though so 6 bits. */
  int rf_base;         /* = 60 */

  int nint_r_g_bits;   /* number of bits in nint_r_g in [60:130] so 8 bits (0:255)
			  we will subytact 50 though so 7 bits. */
  int rg_base;         /* = 60 */  

  int nhuebin_bits;    /* number of bits in huebin   in [0:16]   so 5 bits. */
  int nclass1_bits;    /* = 6+7+5 = 18 */

  int c1_rf_offset;    /* 60 */
  int c1_rg_offset;    /* 60 */

  int irf_mask;        /* 63  */
  int irg_mask;        /* 127 */

  int ihb_mask;        /* 31  */
  int rf_shift;        /* 12  */

  int rg_shift;        /* 5   */
  int c1_num_sh_bins;  /* 36 each of the single cell p, v, or f,bins and 
                          the 27 non-neutral combinations of p,v, and f */

  int nclass2_bits;    /* = 6 */
  int pvf_mask;

  int p_shift;
  int v_shift;

  int f_shift;
  int p_sh_bits;       /* = 17 */

  int p_sh_length;      /* = 2^(5+6+6) = 2^17 = 131072 */
  int p_sh_rf_base;    /* = 70 */

  int p_sh_rf_bits;    /* = 5 [0-30] */
  int p_sh_rf_mask;    /* 31 */

  int p_sh_rf_shift;   /* = 12  */
  int p_sh_rg_base;    /* = 89 */

  int p_sh_rg_bits;    /* = 6 [0:51] */
  int p_sh_rg_mask;    /* = 63 */

  int p_sh_rg_shift;    /* = 6 */
  int p_sh_rcsh1_base;  /* = -12 */

  int p_sh_rcsh1_bits;  /* = 6 [0:35] */
  int p_sh_rcsh1_mask;	/* = 63 */

  int v_sh_bits;        /* = 12 */
  int v_sh_length;      /* = 2^(6+6) = 2^12 = 4096 */
 
  int v_sh_rg_base;     /* = 100 */
  int v_sh_rg_bits;     /* = 6, [0-40] */

  int v_sh_rg_mask;     /* = 63 */
  int v_sh_rg_shift;    /* = 6 */

  int v_sh_rcsh1_base;  /* = 0 */
  int v_sh_rcsh1_bits;  /* = 6 [0:50] */

  int v_sh_rcsh1_mask;  /* = 63 */
  int f_sh_bits;        /* = 8 */

  int f_sh_length;      /* = 2^*(4+4) 256 */
  int f_sh_rf_base;        /* = 85 */

  int f_sh_rf_bits;        /* = 4 [0:15] */
  int f_sh_rf_mask;        /* = 15 */

  int f_sh_rf_shift;       /* = 4 */
  int f_sh_rfh1_base;      /* = 85 */

  int f_sh_rfh1_bits;      /* = 4 [0:15] */
  int f_sh_rfh1_mask;      /* = 15 */

  int number_counts; /* see alloc3 */
  int neutral_counts_offset; /* index of the first of the four neutral group counts in all_counts */

  int number_results_offset; /* index of the number_results count in the all_counts vector. */
  int ix_delta_sum_offset;

  int iy_delta_sum_offset;
  int num_in_tca_offset;

  int num_labels;
  int length_c1_sh; /* = c1_num_sh_bins * length_class1_group */

  int c1_sh_opt;    /* 0 for indented class1 subhistogram
		       1 for very long class1 histogram lines */
  int ipad;
}
;
#endif
