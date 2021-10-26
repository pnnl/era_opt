#include "system_includes.h"
#include "state_struct.h"
#include "matrix_transpose.h"
#include "extract_subset.h"
#include "mtmp.h" /* dgemm customization */
#ifdef TM30
#include "ledfactor.h"
#include "ledinvert.h"
#include "init_huebins.h"
#include "input_mcat02inv.h"
#include "calc_ref_ucs_table.h"
#include "print_cmf.h"
#include "mtvp.h" /* dgemv customization */
#endif
/*
#include "blas.h"
#include "lapack.h"
*/
#include "prep.h"
int prep(struct state_struct *state) {
  /*
    if subset_size < num_channels.
    Extract the subset of columns from the channels matrix replacing it
    with the channel_subset 
    Form the transposes of the channels and source matrices,
    form the channels_avg matrix,
    compute the equations of the edges of the huebins
    and the inside sign indicator for each edge,  (t_a_b_c_s)
    Set pointers to relevant sub-pieces of the source and
    channels matrices.
    Compute the ref_ucs_table of ucs coordinates (japos,mapos,aaposm, baposm, href)
    for each sample at each reference t_t temperature.

    Called by: era_opt
    Calls:     matrix_transpose, extract_subset, init_huebins, mtmp, ledfactor, ledinvert, 
               input_mcat02inv, calc_ref_ucs_table, mtvp, print_cmf,
               exp, log, cos, sin, atan

    Fields of state used:
      *source_matrix,
      *channel_matrix,
      *channel_subset,
      *s026_matrix,
      *huebins
      *subset_indices,
      source_matrix_t,
      channel_matrix_t,
      source_rows,
      source_columns,
      channels_rows,
      channels_columns,
      subset_size,
      s026_columns,
      sbar_t_channels,
      ybar_t_channels,
      ybar_10_t_channels,
      e_t_channels,
      e_t,
      sbar_10_t_channels,
      rx_sbar_10,
      rx_sbar_10_t_channels,
      num_channels,
      num_radiator_temps,
      first_channel_index
      last_channel_index,
      wavelength_index
      xbar_index
      ybar_index
      zbar_index,
      xbar_10_index,
      ybar_10_index,
      zbar_10_index,
      s_0_index,
      s_1_index,
      s_2_index,
      hue_bins,
      mcat02,
      mcat02_lu,
      mcat02_ipiv,
      mhpe,
      tt,
      wavelength,
      table_planck,
      read_mcat02inv,
      sref_mode,

    Fields of state set:
      channels,
      channel_subset,
      t_a,
      t_b,
      t_c,
      t_s,
      wavelength,
      xbar,
      ybar,
      zbar,
      xbar_10,
      ybar_10,
      zbar_10,
      s_0,
      s_1,
      s_2,
      *source_matrix_t,
      *channel_matrix_t
      *s026,
      *huebins,
      *sbar_t_channels,
      *ybar_t_channels,
      *e_t,
      *e_t_channels,
      *sz,
      *sz_t_channels,
      *sbar_10_t_channels,
      *rx_sbar_10,
      *rx_sbar_10_t_channels,
      *t_a_b_c_s,
      *mcat02
      *mcat02_lu
      *mcat02_ipiv
      *mhpe
      *tt,
      *wavelength,
      *table_planck,
      *chroma_shift_rr,
      fl,
      fl_over_100
      fourth_root_fl;
      zz,
      nbb,
      ncc
      pi,
      degrees_per_radian,
      radians_per_degree,
      c,
      f,
      sr,
      nc,
      chroma_c_scale,
      ref_ucs_table
  */
  double *source_matrix;
  double *channel_matrix;
  double *channel_subset;
  double *s026_matrix;
  /*
  double *quadrangles;
  */
  double *huebins;
  double *samples;
  double *sbar_t_channels;
  double *ybar_t_channels;
  double *e_t_channels;
  double *sz_t_channels;
  double *e_t;
  double *sz;
  double *sbar_10_t_channels;
  double *ybar_10_t_channels;
  double *rx_sbar_10;
  double *rx_sbar_10_t_channels;
  double *channels;
  double *s026;
  double *wavelength;
  double *xbar;
  double *ybar;
  double *zbar;
  double *sbar;
  double *xbar_10;
  double *ybar_10;
  double *zbar_10;
  double *sbar_10;
  double *s_0;
  double *s_1;
  double *s_2;
  double *source_matrix_t;
  double *channel_matrix_t;
  double *tt;
  double *table_planck;
  double *mu_wavelength;
  double *tplanck_num;
  double *c2_by_wavelength;
  double *tplanck;
  double *rx_sbar_10_col;
  double *rx_col;
  double xyz_planck[3];
  /*
  double *q_a_b_c;
  */
  double *t_a_b_c_s;
  /*
  double *q_a;
  double *q_b;
  double *q_c;
  double *q_s;
  */
  double *mcat02;
  double *mcat02_lu;
  double *mcat02_i;
  double mcat02_i_t[9];
  double *mhpe;
  double *mhpe_mcat02_i;
  double *table_pos;
  double *chroma_shift_rr;
  double *stc_p;
  int    *mcat02_ipiv;
  int    *subset_indices;
  double alpha;
  double beta;
  double fl;
  double zz;
  double nbb;
  double ncb;
  double la;
  double five_la;
  double cube_rt_five_la;
  double k;
  double k_to_4;
  double one_m_k_to_4;
  double yb;
  double yw;
  double n;
  double fifth_rt_recip_n;
  double c;
  double chroma_c_scale;
  double f;
  double sr;
  double nc;
  double lambda;
  double mu_lambda;
  double mu_lambda_5;
  double recip_tt;
  double max_tplanck;
  double recip_max_tplanck;
  double tplanck_kk;
  double ui_vi_denom;
  double ui;
  double vi;
  double c1;
  double c2;
  double pi_by_num_huebins;
  double cangle;
  double cscale;
  double fl_over_100;
  double fourth_root_fl;
  double half_pi;
  double pi;
  double two_pi;
  double degrees_per_radian;
  double radians_per_degree;
  double x;

  double twentieth;
  double eleventh;
  double ninth;
  double twelve_elevenths;
  double sin_2;
  double cos_2;
  double t_const;
  double degree_adaptation;
  double hr_thresh;
  
  int char_len;
  int ithree;

  int source_rows;
  int source_columns;

  int channels_rows;
  int channels_columns;

  int num_channels;
  int num_radiator_temps;

  int first_channel_index;
  int last_channel_index;

  /*
  int wavelength_index;
  */
  int xbar_index;
  int xbar_10_index;
  /*
  int ybar_index;
  int zbar_index;
  */

  int s_0_index;
  int rad_inc;

  /*
  int ybar_10_index;
  int zbar_10_index;
  */

  /*
  int s_1_index;
  int s_2_index;
  */

  int i;
  int j;

  int three_num_samples;
  int num_samples;

  int info;
  int kk;

  int inc1;
  int two_num_huebins;

  int success;
  int num_huebins;

  int s026_columns;
  int sz_index;
  
  int ipos;
  int jpos;

  int kpos;
  int subset_size;

  int read_mcat02inv;
  int sref_mode;

  int select_only;

  char t_char;
  char n_char;

  FILE *lfp;
  FILE *efp;
  
  /*
  c1=3.7417749/1e16;
  */
  c1 =3.741832/1e16;
  state->c1           = c1;
  char_len = 1;
  inc1     = 1;
  ithree = 3;
  t_char = 'T';
  n_char = 'N';
  lfp                 = state->lfp;
  channel_matrix      = state->channel_matrix;
  channel_matrix_t    = state->channel_matrix_t;
  channel_subset      = state->channel_subset;
  select_only         = state->select_only;
  subset_size         = state->subset_size;
  subset_indices      = state->subset_indices;
  source_matrix       = state->source_matrix;
  source_matrix_t     = state->source_matrix_t;
  sbar_t_channels     = state->sbar_t_channels;
  source_rows         = state->source_rows;
  source_columns      = state->source_columns;
  channels_rows       = state->channels_rows;
  channels_columns    = state->channels_columns;
  first_channel_index = state->first_channel_index;
  last_channel_index  = state->last_channel_index;
  num_channels        = state->num_channels;
  num_samples         = state->num_samples;
  num_huebins         = state->num_huebins;
  sz_index            = state->sz_index;
  /*
    wavelength_index    = state->wavelength_index;
  */
  xbar_index          = state->xbar_index;
  /*
    ybar_index          = state->ybar_index;
    zbar_index          = state->zbar_index;
  */
  xbar_10_index       = state->xbar_10_index;
  /*
    ybar_10_index       = state->ybar_10_index;
    zbar_10_index       = state->zbar_10_index;
  */
  s_0_index           = state->s_0_index;
  /*
    s_1_index           = state->s_1_index;
    s_2_index           = state->s_2_index;
  */
  /*
    wavelength is needed in gen_combos to generate the subset channels file.
  */
  wavelength          = state->wavelength;

  if (select_only == 0) {
#ifdef TM30
    s026_matrix         = state->s026_matrix;
    /*
      quadrangles         = state->quadrangles;
    */
    huebins             = state->huebins;
    sbar_10_t_channels  = state->sbar_10_t_channels;
    rx_sbar_10          = state->rx_sbar_10;
    rx_sbar_10_t_channels = state->rx_sbar_10_t_channels;
    s026_matrix         = state->s026_matrix;
    s026                = state->s026;
    sz_t_channels       = state->sz_t_channels;
    /*
      q_a_b_c             = state->q_a_b_c;
      q_s                 = state->q_s;
    */
    tt                  = state->tt;
    table_planck        = state->table_planck;
    t_a_b_c_s           = state->t_a_b_c_s;
    s026_columns        = state->s026_columns;
    mcat02              = state->mcat02;
    mcat02_lu           = state->mcat02_lu;
    mcat02_i            = state->mcat02_i;
    mcat02_ipiv         = state->mcat02_ipiv;
    mhpe                = state->mhpe;
    mhpe_mcat02_i       = state->mhpe_mcat02_i;
    mu_wavelength       = state->mu_wavelength;
    tplanck_num         = state->tplanck_num;
    c2_by_wavelength    = state->c2_by_wavelength;
    tplanck             = state->tplanck;
    chroma_shift_rr     = state->chroma_shift_rr;
    c1                  = state->c1;
    num_radiator_temps  = state->num_radiator_temps;
    ybar_t_channels     = state->ybar_t_channels;
    ybar_10_t_channels  = state->ybar_10_t_channels;
    e_t                 = state->e_t;
    e_t_channels        = state->e_t_channels;
    read_mcat02inv      = state->read_mcat02inv;
    sref_mode           = state->sref_mode;
#endif
  }
  /*
    Tanspose the channels  matrix
  */
  matrix_transpose(channels_rows,channels_columns,
		   channel_matrix,channel_matrix_t);
  num_channels = last_channel_index - first_channel_index + 1;
  /*
    Caution address arithmetic.
  */
  channels = channel_matrix_t + ((first_channel_index - 1) * channels_rows);
  state->num_channels = num_channels;
  state->channels     = channels;
  state->all_channels = channels;
  /*
    If subset_size < num_channels extract the subset based on the
    subset_indices vector.
  */
  if (subset_size < num_channels) {
    extract_subset(channels_rows,subset_size,subset_indices,channels,channel_subset);
    /*
      reset channels and num_channels to the channel_subset and subset_size
      respectively, the original channels pointer is in state->all_channels,
      and the original num_channels is in state->full_num_channels - see
      read_params.c
    */
    channels            = channel_subset;
    num_channels        = subset_size;
    state->num_channels = subset_size;
    state->channels     = channel_subset;
  }
  /*
    Tanspose the source matrix
  */
  matrix_transpose(source_rows,source_columns,source_matrix,source_matrix_t);
  /*
    Set pointers to the unit stride columns of the source matrix in
    its transpose.
  */
  samples     = source_matrix_t;
  /* 
     Caution the following 6 statments do address arithmetic
  */
  xbar        = source_matrix_t + ((xbar_index - 1) * source_rows);
  ybar        = xbar + source_rows;
  zbar        = ybar + source_rows;
  xbar_10     = source_matrix_t + ((xbar_10_index - 1) * source_rows);
  ybar_10     = xbar_10 + source_rows;
  zbar_10     = ybar_10 + source_rows;
  /*
    wavelength = source_matrix_t + ((wavelength_index - 1) * source_rows);
  */
  sbar        = xbar; 
  sbar_10     = xbar_10;
  state->sbar        = sbar;
  state->xbar        = xbar;
  state->ybar        = ybar;
  state->zbar        = zbar;
  state->sbar_10     = sbar_10;
  state->xbar_10     = xbar_10;
  state->ybar_10     = ybar_10;
  state->zbar_10     = zbar_10;
  /*
    sbar_t_channels is a 3 x num_channels array stored in column major order
    sbar is source_rows x 3 stored in column major order.
    channels is source_rows x num_channels stored in column major order.
    sbar_t_channels = sbar'*channels and is  3 x num channels stored in 
    column major order.
  */
  success = mtmp(source_rows,ithree,sbar,
		 source_rows,num_channels,channels,sbar_t_channels,lfp);

  half_pi = 2.0*atan(1.0);
  pi = 2.0*half_pi;
  /*
    pi = 4.0*atan(1.0);
  */
  two_pi = 2.0 * pi;
  degrees_per_radian = 180.0/pi;
  radians_per_degree = pi/180.0;

  state->half_pi            = half_pi;
  state->pi                 = pi;
  state->two_pi             = two_pi;
  state->radians_per_degree = radians_per_degree;
  state->degrees_per_radian = degrees_per_radian;
  /*
    Set the wavelengths vector.
  */
  j = 0;
  for (i=380;i<=780;i+=1) {
    lambda = (double)i;
    wavelength[j] = lambda;
    j += 1;
  }

  if (select_only == 0) {
#ifdef TM30
#ifdef DBG
    print_cmf(source_rows,xbar,ybar,zbar,xbar_10,ybar_10,zbar_10,lfp);
#endif
    /*
      Caution the following 3 statments do address arighmetic.
    */
    s_0         = source_matrix_t + ((s_0_index - 1) * source_rows);
    s_1         = s_0 + source_rows;
    s_2         = s_1 + source_rows;
    state->s_0         = s_0;
    state->s_1         = s_1;
    state->s_2         = s_2;
    /* 
      sbar is a source_rows x 3 matrix stored in column major order 

      Transpose the s026_matrix
    */
    matrix_transpose(source_rows,s026_columns,s026_matrix,s026);
    sz = s026 + ((sz_index-1)*source_rows); /* Caution address arithmetic */
    state->sz = sz;
    /*
      channels is channels_rows x num_channels 
      stored in column major order.
      we want to compute sbar_t_channels = s_avg^T * channels
    */
    /*
      alpha = 1.0;
      beta  = 0.0;
      dgemm_(&t_char,&n_char,&ithree,&num_channels,&source_rows,
             &alpha,sbar,&source_rows,
             channels,&source_rows,&beta,
	     sbar_t_channels,&ithree,char_len,char_len);
    */
    /*
      Extract the ybar_t_channels vector from sbar_t_channels.
    */
    stc_p = &sbar_t_channels[1];
    for (i=0;i<num_channels;i++) {
      ybar_t_channels[i] = *stc_p;
      stc_p += 3; /* Caution address arithmetic */
    }
    for (i=0;i<source_rows;i++) {
      e_t[i] = 1.0;
    }

    mtvp(source_rows,num_channels,channels,source_rows,e_t,e_t_channels);
    /*
      alpha = 1.0;
      beta  = 0.0;
      dgemv_(&t_char,&source_rows,&num_channels,&alpha,channels,&source_rows,
      e_t,&inc1,&beta,e_t_channels,&inc1,char_len);
    */
    success = mtmp(source_rows,ithree,sbar_10,
		   source_rows,num_channels,channels,sbar_10_t_channels,lfp);
    /*
      alpha = 1.0;
      beta  = 0.0;
      dgemm_(&t_char,&n_char,&ithree,&num_channels,&source_rows,&alpha,sbar_10,
             &source_rows,channels,&source_rows,&beta,
             sbar_10_t_channels,&ithree,char_len,char_len);
    */
    /*
      Extract ybar_10_t_channels from sbar_10_t_channels
    */
    stc_p = &sbar_10_t_channels[1];
    for (i=0;i<num_channels;i++) {
      ybar_10_t_channels[i] = *stc_p;
      stc_p += 3; /* Caution address arithmetic */
    }

    mtvp(source_rows,num_channels,channels,source_rows,sz,sz_t_channels);
    /*
      alpha = 1.0;
      beta = 0.0;
      dgemv_(&t_char,&source_rows,&num_channels,&alpha,channels,&source_rows,
   	     sz,&inc1,&beta,sz_t_channels,&inc1,char_len);
    */


    /*
      sbar_10_t_channels is a 3 x num_channels array stored in column major order
      (from the fortran dgemm_ call)
    */
    /*
      Form the rx_sbar_10 matrix which is the elementwise product of each sample 
      column with the xbar_10, ybar_10, zbar_10 vectors
      in succession the result is a source_rows x 3*num_samples matrix stored
      in column-major order.
    */
    rx_sbar_10_col = rx_sbar_10;
    rx_col         = samples;
    for (i=0;i<num_samples;i++) {
      for (j=0;j<source_rows;j++) {
	rx_sbar_10_col[j] = rx_col[j] * xbar_10[j];
      }
      rx_sbar_10_col += source_rows; // Caution address arithmetic 
      for (j=0;j<source_rows;j++) {
	rx_sbar_10_col[j] = rx_col[j] * ybar_10[j];
      }
      rx_sbar_10_col += source_rows; // Caution address arithmetic 
      for (j=0;j<source_rows;j++) {
	rx_sbar_10_col[j] = rx_col[j] * zbar_10[j];
      }
      rx_sbar_10_col += source_rows; // Caution address arithmetic 
      rx_col += source_rows; // Caution address arithmetic 
    } // end for (i... *)
    /*
      form the rx_sbar_10_t_channels matrix which is 
      3*num_samples x num_channels and is rx_sbar_10' * channels
      stored in column major order.
    */
    three_num_samples = 3*num_samples;
    success = mtmp(source_rows,three_num_samples,rx_sbar_10,
		   source_rows,num_channels,channels,
		   rx_sbar_10_t_channels,lfp);
    /*
    alpha = 1.0;
    beta  = 0.0;
    dgemm_(&t_char,&n_char,&three_num_samples,&num_channels,&source_rows,
           &alpha,rx_sbar_10,&source_rows,channels,&source_rows,&beta,
	   rx_sbar_10_t_channels,&three_num_samples,char_len,char_len);
    */
    /*
      Determine q_a_b_c, and q_s, the coefficients and side signs of the
      edges of the quadrilaterals with coordinates in quadrangles
      success = init_ansi_check(num_quadrangles,quadrangles,q_a_b_c,q_s);
    */
    /*
      Determine the t_a_b_c_s, the coefficients and side signs of the
      edges of the triangles with coorodinates in triangles.
    */
    success = init_huebins(num_huebins,huebins,t_a_b_c_s,lfp);
    /*
      Initialize the mcat02, mcat02_lu and mcat02_ipiv arrays.
      MCAT02=[0.7328 0.4296 -0.1624;-0.7036 1.6975 0.0061;0.0030 0.0136 0.9834];
    */
    if (success) {
      mcat02[0] =  0.7328;
      mcat02[1] = -0.7036;
      mcat02[2] =  0.0030;
      mcat02[3] =  0.4296;
      mcat02[4] =  1.6975;
      mcat02[5] =  0.0136;
      mcat02[6] = -0.1624;
      mcat02[7] =  0.0061;
      mcat02[8] =  0.9834;
      if (read_mcat02inv == 0) {
	for (i=0;i<9;i++) {
	  mcat02_lu[i] = mcat02[i];
	}
	ledfactor(mcat02_lu);
	ledinvert(mcat02_lu,mcat02_i);
      }
#ifdef DBG
      if (lfp) {
	fprintf(lfp,"mcat02^-1\n");
	fprintf(lfp,"%17.14f %17.14f %17.14f\n",
		mcat02_i[0],mcat02_i[3],mcat02_i[6]);
	fprintf(lfp,"%17.14f %17.14f %17.14f\n",
		mcat02_i[1],mcat02_i[4],mcat02_i[7]);
	fprintf(lfp,"%17.14f %17.14f %17.14f\n",
		mcat02_i[2],mcat02_i[5],mcat02_i[8]);
	fflush(lfp);
      }
#endif
      if (read_mcat02inv != 0) {
	success = input_mcat02inv(state);
      }
      /*
	ithree = 3;
	dgetrf_(&ithree,&ithree,mcat02_lu,&ithree,mcat02_ipiv,&info);
	if (info != 0) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"prep: mcat02 matrix was singular\n");
        }
      }
      */
    }
    if (success) {
      /*
	MHPE=[0.38971 0.68898 -0.07868;-0.22981 1.18340 0.04641;0 0 1];
      */
      mhpe[0] =  0.38971;
      mhpe[1] = -0.22981;
      mhpe[2] =  0.0;
      mhpe[3] =  0.68898;
      mhpe[4] =  1.18340;
      mhpe[5] =  0.0;
      mhpe[6] = -0.07868;
      mhpe[7] =  0.04641;
      mhpe[8] =  1.0;
      /*
	Now we form mhpe * mcat02^-1
	This is a 3x3 matrix matrix multiply.
	as in test_ledla.c
      */
      ipos = 0;
      kpos = 0;
      for (j=0;j<3;j++) {
	/*
	  compute column j of mhpe_mcat02_i. mhpe_mcat02_i[ipos:ipos+2]
	  loop over elements in column j of mcat02_i;
	*/
	jpos = 0;
	for (kk=0;kk<3;kk++) {
	  x = mcat02_i[kpos+kk];
	  for (i=0;i<3;i++) {
	    mhpe_mcat02_i[ipos+i] += x * mhpe[jpos+i];
	  }
	  /*
	    Move to the next column of mcat02.
	    and next element of mcat02_i[k]
	  */
	  jpos += 3;
	}
	kpos += 3;
	ipos += 3;
      }
    }
    /*
      Planckian radiation constants
    */
    c2=0.014388;  
    /*
      Set constants used by cam02ucs
    */
    la = 100.0;
    five_la = 5.0*la;
    k=1.0/(five_la+1.0);
    cube_rt_five_la = exp(log(five_la)/3.0);
    k_to_4 = k * k;
    k_to_4 = k_to_4 * k_to_4;
    one_m_k_to_4 = 1.0 - k_to_4;
    fl = (.2 * k_to_4*five_la ) + 
    	(.1 * one_m_k_to_4 * one_m_k_to_4 * cube_rt_five_la);
    fl_over_100 = fl*.01;
    fourth_root_fl = sqrt(fl);
    fourth_root_fl = sqrt(fourth_root_fl);

    yb = 20;
    yw = 100;
    n  = yb/yw;
    zz = 1.48 + sqrt(n);
    fifth_rt_recip_n = exp(log(1.0/n)/5.0);
    nbb = 0.725 * fifth_rt_recip_n;
    ncb = nbb;
    sr = 1.0;
    f  = 1.0;
    c  = 0.69;
    nc = 1.0;
    chroma_c_scale = exp(log(1.64 - exp(log(.29)*n))*.73);
    /*
      Set constants used in the new cam02ucs routine.
    */
    twentieth = 0.05;
    eleventh  = 1.0/11.0;
    ninth     = 1.0/9.0;
    twelve_elevenths = 12.0 * eleventh;
    cos_2      = cos(2.0);
    sin_2      = sin(2.0);
    hr_thresh  = 20.14 * radians_per_degree;
    t_const    = (50000/13.0)*nc*ncb;
    degree_adaptation = f * (1.0-(1.0/3.6)*exp((-la-42)/92));
    state->n                  = n;
    state->yb                 = yb;
    state->fl                 = fl;
    state->fl_over_100        = fl_over_100;
    state->fourth_root_fl     = fourth_root_fl;
    state->zz                 = zz;
    state->nbb                = nbb;
    state->ncb                = ncb;
    state->c                  = c;
    state->c_zz               = c * zz;
    state->f                  = f;
    state->sr                 = sr;
    state->nc                 = nc;
    state->chroma_c_scale     = chroma_c_scale;
    state->c1                 = c1;
    state->c2                 = c2;

    state->twentieth = twentieth;
    state->twelve_elevenths = twelve_elevenths;
    state->eleventh  = eleventh;
    state->ninth     = ninth;
    state->cos_2     = cos_2;
    state->sin_2     = sin_2;
    state->hr_thresh = hr_thresh;
    state->t_const   = t_const;
    state->degree_adaptation = degree_adaptation;
    /*
      Set the tt, wavelengths, ui and vi vectors for planck_rad;
    */
    i = 1000;
    rad_inc = 20;
    for (j=0;j<num_radiator_temps;j++) {
      tt[j] = (double)i;
      i += rad_inc;
    }
    /*
    for (i=1000;i<=20000;i+=20) {
      tt[j] = (double)i;
      j += 1;
    }
    */
    j = 0;
    c2=0.014388;
    /*
      Should change this to take the lambda's from the first column of the
      channels matrix. Best would be to compare for equality to the first
      column of the source_t matrix, and if not equal, interpolate/extrapolate
      the channels data to match the source matrix wavelengths.
      Currently if source_rows != channels_rows we print an error message
      and quit.
    for (i=0;i<channels_rows;i++) {
      lambda = channel_matrix_t[i];
    */
    for (i=380;i<=780;i+=1) {
      lambda = (double)i;
      mu_lambda = lambda * 0.000000001;
      mu_lambda_5 = mu_lambda * mu_lambda;
      mu_lambda_5 = mu_lambda_5 * mu_lambda_5;
      mu_lambda_5 = mu_lambda_5 * mu_lambda;
      mu_wavelength[j] = mu_lambda;
      tplanck_num[j] = c1/mu_lambda_5;
      c2_by_wavelength[j] = c2/mu_lambda;
      j += 1;
    }
    table_pos = table_planck;
    for (i=0;i<num_radiator_temps;i++) {
      *table_pos = tt[i];
      table_pos += 1; /* Caution address arithmetic */
      recip_tt = 1.0/tt[i];
      
      max_tplanck = tplanck_num[0]/(exp(c2_by_wavelength[0]*recip_tt)-1.0);
      tplanck[0] = max_tplanck;
      for (kk=1;kk < source_rows;kk++) {
	tplanck_kk  = tplanck_num[kk]/(exp(c2_by_wavelength[kk]*recip_tt)-1.0);
	tplanck[kk] = tplanck_kk;
	if (tplanck_kk > max_tplanck) {
	  max_tplanck = tplanck_kk;
	}
      }
      recip_max_tplanck = 1.0/max_tplanck;
      for (kk = 0;kk< source_rows;kk++) {
	tplanck[kk] *= recip_max_tplanck;
      }
      alpha = 1.0;
      beta  = 0.0;
      
      mtvp(source_rows,ithree,sbar,source_rows,tplanck,xyz_planck);
      /*
      dgemv_(&t_char,&source_rows,&ithree,&alpha,sbar,&source_rows,
	     tplanck,&inc1, &beta,xyz_planck,&inc1,char_len);
      */
      ui_vi_denom = 1.0/(xyz_planck[0] + (15.0*xyz_planck[1]) + 
			 (3.0 *xyz_planck[2]));
      ui = 4.0 * xyz_planck[0] * ui_vi_denom;
      vi = 6.0 * xyz_planck[1] * ui_vi_denom;
      *table_pos = ui;
      table_pos += 1; /* Caution address arithmetic */
      *table_pos = vi;
      table_pos += 2; /* Caution address arithmetic, skip di_sq spot
		         as that is filled in in planck_rad routine */
    } /* end for i */
    if (success) {
      /*
	fill the chroma_shift_rr matrix. Have this be a 
	num_huebins x 2 matrix stored in row major order
	with the row i being 100.0 * cos(2i+1*pi/16), 100.0 * sin(2i+1*pi/16)
	for i= 0..num_huebins-1
	Note this matrix is very much tied to the specific huebins
	being used, if you change those you will want to change
	these perhaps they should be read in
      */
      pi_by_num_huebins = pi/((double)num_huebins);
      j = 0;
      for (i=0;i<num_huebins;i++) {
	cangle = ((2.0*i)+1.0)*pi_by_num_huebins;
	chroma_shift_rr[j] = cos(cangle);
	chroma_shift_rr[j+1] = sin(cangle);
	j += 2;
      }
      cscale = 100.0;
      two_num_huebins = num_huebins + num_huebins;
      for (i=0;i<two_num_huebins;i++) {
	chroma_shift_rr[i] = chroma_shift_rr[i]*cscale;
      }
      /*
	dscal_(&two_num_huebins,&cscale,chroma_shift_rr,&inc1);
      */
      /*
	Variables for tracking t_t ranged from calc_cct and max rf value.
      */
      state->cct_t_t_min = 4000.0;
      state->cct_t_t_max = 3000.0;
      state->rf_max  = 0.0;
    }
    if (success) {
      /*
	Set the ref_ucs_table if sref_mod > 0.
      */
      if (sref_mode > 0) {
	success = calc_ref_ucs_table(state);
      }
    }
#endif /* TM30 */
  } /* end if (select_only == 0) */
  return(success);
}
