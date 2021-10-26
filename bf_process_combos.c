#include "system_includes.h"
#include "state_struct.h"
#include "unpack_isteps.h"
#include "count_nzc.h"
#include "process_combo.h"
/*
#include "unique_ratio.h"
#include "get_steps.h"
#include "ies_tm_30_18.h"
#include "output_tm_30_results.h"
#include "output_tm_30_debug.h"
*/
#include "bf_process_combos.h"
void bf_process_combos(struct state_struct *state,
		       struct loop_state_struct *loop_state, 
		       int myid, int nproc) {
  /*
    Triply nested set of pnnl loops to test all combinations.
    Called by: era_opt
    Calls:     unpack_isteps, count_nzc, process_combo
  */
  double *fracts;
  double *tm_30_results;
  double *steps;
  double *test_radii_sq;
  double *ybar_t_channels;
  double *e_t_channels;
  double *sz_t_channels;
  double *z_dom_data;    /* state field */
  double *y_dom_data;   /* state field */
  double *x_dom_data;     /* state field */

  double *my_z_dom_data;
  double *my_y_dom_data;
  double *my_x_dom_data;
  double cct_vals[4];

  int64_t *iy_dom_bins; /* state field */
  int64_t *ix_dom_bins;   /* state field */

  int *istep_index;
  int *mm_istep_index;
  int *imap_sort;
  int *unique_prime_factors;
  int *upf_index;
  int *iz_dom_steps;
  int *iy_dom_steps;
  int *ix_dom_steps;

  double z;
  double y;
  double x;
  double uprime;
  double vprime;
  double denom;
  double rdenom;

  double udist;
  double vdist;
  double u_target;
  double v_target;
  double sqdist;
  double rf_thresh;
  double ybar_dot_stest;
  double stest_sum;
  double ler;
  double z_dom_x;
  double z_dom_y;
  double z_dom_z;
  double y_dom_x;
  double y_dom_y;
  double y_dom_z;
  double x_dom_x;
  double x_dom_y;
  double x_dom_z;
  double teal_x;
  double teal_y;
  double teal_z;
  double sum_m;
  double sum_p;
  double m_over_p;
  double dist;

  double watchme;
  double progress;
  double pdone;
  double precip_izc;

  int64_t izc;

  int64_t iy_start;
  int64_t iy_end;
  int64_t ix_start;
  int64_t ix_end;

  int64_t iz;
  int64_t iy;
  int64_t ix;

  int64_t iy_delta_sum;
  int64_t ix_delta_sum;
  int64_t lix_delta_sum;
  int64_t iy_delta;
  int64_t ix_delta;
  int64_t num_in_tca;
  int64_t max_iy_delta;
  int64_t max_ix_delta;

  int num_channels;
  int full_output;

  int nz;
  int ny;

  int i;
  int nx;

  int nzc;
  int extra_log_info;

  int nz_dom_fields;     /* state field */
  int iz_dom_tsv_offset; /* state field */

  int iz_dom_delta; 
  int iz_dom_sbpi;       /* state field */

  int ny_dom_fields;    /* state field */
  int iy_dom_tsv_offset; /* state field */

  int iy_dom_sbpi;       /* state field */
  int nx_dom_fields;       /* state field */

  int ix_dom_tsv_offset;   /* state field */
  int ix_dom_sbpi;         /* state field */

  int ny_dom_bins;
  int nx_dom_bins;

  int ivb_low;
  int ivb_high;

  int iub_low;
  int iub_high;

  int nzzc;
  int nzyc;
  
  int nztc;
  int nzxc;

  int unique;
  int max_index;


  FILE *lfp;
  FILE *ofp;

  istep_index  	      = loop_state->istep_index;
  mm_istep_index      = loop_state->mm_istep_index;
  imap_sort           = loop_state->imap_sort;
  unique_prime_factors = loop_state->unique_prime_factors;
  upf_index           = loop_state->upf_index;
  max_index           = loop_state->max_index;
  num_channels 	      = loop_state->num_channels;
  fracts              = loop_state->curr_xyz;
  nz                  = loop_state->nz;
  ny                  = loop_state->ny;
  nx                  = loop_state->nx;
  lfp                 = loop_state->lfp;
  fracts              = loop_state->fracts;
  tm_30_results       = state->tm_30_results;
  rf_thresh           = state->rf_thresh;
  ofp                 = state->ofp;
  extra_log_info      = state->extra_log_info;
  test_radii_sq       = state->test_radii_sq;
  ybar_t_channels     = state->ybar_t_channels;
  e_t_channels        = state->e_t_channels;
  sz_t_channels       = state->sz_t_channels;
  full_output         = state->full_output;
  watchme             = state->watchme;
  izc = 0;
  u_target = loop_state->u_target;
  v_target = loop_state->v_target;
  /*
  leeway = 0,01;
  mratio = 1.0-leeway;
  pratio = 1.0+leeway;
  xmin_c = xmin_c * mratio;
  ymin_c = ymin_c * mratio;
  xmax_c = xmax_c * pratio;
  ymax_c = ymax_c * pratio;
  */
  
  iz_dom_sbpi        = state->iz_dom_sbpi;
  iz_dom_tsv_offset  = state->iz_dom_tsv_offset;
  nz_dom_fields      = state->nz_dom_fields;
  izc               = state->izc;
  if (izc != 0) {
    precip_izc = 100.0/((double)izc);
  } else {
    precip_izc = 0.0;
  }
  z_dom_data         = state->z_dom_data;
                   
  iy_dom_sbpi       = state->iy_dom_sbpi;
  iy_dom_tsv_offset = state->iy_dom_tsv_offset;
  ny_dom_fields     = state->ny_dom_fields;
  ny_dom_bins       = state->ny_dom_bins;
  y_dom_data        = state->y_dom_data;
  iy_dom_bins       = state->iy_dom_bins;
                   
  ix_dom_sbpi         = state->ix_dom_sbpi;
  ix_dom_tsv_offset   = state->ix_dom_tsv_offset;
  nx_dom_fields       = state->nx_dom_fields;
  nx_dom_bins         = state->nx_dom_bins;
  x_dom_data          = state->x_dom_data;
  ix_dom_bins         = state->ix_dom_bins;

  iy_delta_sum = (int64_t)0;
  ix_delta_sum = (int64_t)0;
  num_in_tca   = (int64_t)0;
  max_iy_delta = (int64_t)0;
  max_ix_delta = (int64_t)0;
  iz_dom_delta = nz_dom_fields * nproc;
  my_z_dom_data = z_dom_data + (myid * nz_dom_fields); /* Caution, Address arithmetic */
  /*
    Round robin scheduling of tasks.
  */
  progress = 0.0;
  for (iz = myid; iz < izc; iz += nproc) {
    iz_dom_steps = (int*)my_z_dom_data;
    unpack_isteps(nz,iz_dom_sbpi,istep_index,iz_dom_steps);
    nzzc   = count_nzc(nz,istep_index);
    z_dom_x = my_z_dom_data[iz_dom_tsv_offset];
    z_dom_y = my_z_dom_data[iz_dom_tsv_offset+1];
    z_dom_z = my_z_dom_data[iz_dom_tsv_offset+2];
    my_z_dom_data += iz_dom_delta; /* Caution, Address arithmetic. */
    ivb_low = 0;
    ivb_high = ny_dom_bins-1;
    iy_start = iy_dom_bins[ivb_low];
    iy_end   = iy_dom_bins[ivb_high];
    iy_delta = iy_end - iy_start;
    if (iy_delta > max_iy_delta) {
      max_iy_delta = iy_delta;
    }
    iy_delta_sum += iy_delta;
    /*
      Subtract off the existing y.
    */
    /*
      So we want to look at y combos with y between ymin and ymax.
    iymin_bin = (ymin/combo_bin_width);
    iymax_bin = ((ymax/combo_bin_width) + .5);
    iy_start  = iy_dom_bins[iymin_bin];
    iy_end    = iy_dom_bins[iymax_bin];
    if (iy_end >= iyc) {
      iy_end = iyc - 1;
    }
    */
    lix_delta_sum = 0;
    my_y_dom_data = y_dom_data + (iy_start * ny_dom_fields); /* Caution Address arithmetic */
    for (iy = iy_start;iy<=iy_end;iy++) {
      iy_dom_steps = (int*)my_y_dom_data;
      unpack_isteps(ny,iy_dom_sbpi,&istep_index[nz],iy_dom_steps);
      nzyc  = count_nzc(ny,&istep_index[nz]);
      nztc  = nzzc + nzyc;
      y_dom_x = my_y_dom_data[iy_dom_tsv_offset];
      y_dom_y = my_y_dom_data[iy_dom_tsv_offset+1];
      y_dom_z = my_y_dom_data[iy_dom_tsv_offset+2];
      teal_x = z_dom_x + y_dom_x;
      teal_y = z_dom_y + y_dom_y;
      teal_z = z_dom_z + y_dom_z;
      my_y_dom_data += ny_dom_fields; /* Caution Address arithmetic */
      /*
	Now for this blue and green combination we want to examine 
	all of the red combinations 
      */
      iub_low = 0;
      iub_high = nx_dom_bins - 1;
      /*
#define DBG 1
      */
#ifdef DBG
      fprintf(stdout,"iub_low = %d,iub_high=%d\n",iub_low,iub_high);
      fflush(stdout);
#endif
      ix_start = ix_dom_bins[iub_low];
      ix_end   = ix_dom_bins[iub_high];
      ix_delta = ix_end - ix_start + 1;
      lix_delta_sum += ix_delta;
      my_x_dom_data = x_dom_data + (ix_start * nx_dom_fields); /* Cautions Addersss arithmetic */
#ifdef DBG
      fprintf(stdout,"ix_start = %ld, ix_end = %ld\n",ix_start,ix_end);
      fflush(stdout);
#endif
      for (ix = ix_start;ix<=ix_end;ix++) {
	ix_dom_steps = (int*)my_x_dom_data;
	unpack_isteps(nx,ix_dom_sbpi,&istep_index[nz+ny],ix_dom_steps);
	nzxc  = count_nzc(nx,&istep_index[nz+ny]);
	nzc   = nztc + nzxc;
	x_dom_x = my_x_dom_data[ix_dom_tsv_offset];
	x_dom_y = my_x_dom_data[ix_dom_tsv_offset+1];
	x_dom_z = my_x_dom_data[ix_dom_tsv_offset+2];
	x = x_dom_x + teal_x;
	y = x_dom_y + teal_y;
	z = x_dom_z + teal_z;
	my_x_dom_data += nx_dom_fields; /* Caution Address arithmetic */
	denom = x + (15.0*y) + (3.0*z);
	rdenom = 1.0/denom;
	uprime = 4.0*x*rdenom;
	vprime = 9.0*y*rdenom;
	udist = uprime - u_target;
	vdist = vprime - v_target;
	sqdist = (udist * udist) + (vdist * vdist);

	if (sqdist < test_radii_sq[nzc]) {
	  process_combo(state,loop_state,nzc,&num_in_tca,
			istep_index,x,y,z,sqdist);
	} /* end if in target circle */
      } /* end for ix */
    } /* end for iy */
    pdone = iz*precip_izc;
    if (pdone >= progress) {
      if (lfp) {
	fprintf(lfp,"iz = %ld completed, %% done = %5.1f\n",
		iz,pdone);
	fflush(lfp);
      }
      progress = progress+watchme;
    }
    if (lix_delta_sum > max_ix_delta) {
      max_ix_delta = lix_delta_sum;
    }
    ix_delta_sum += lix_delta_sum;
  } /* end for iz */
  if (lfp) {
    fprintf(lfp,"100%% complete\n");
    fprintf(lfp,"iy_delta_sum = %ld\n",iy_delta_sum);
    fprintf(lfp,"ix_delta_sum = %ld\n",ix_delta_sum);
    fprintf(lfp,"num in tca   = %ld\n",num_in_tca);
    fprintf(lfp,"max_iy_delta = %ld\n",max_iy_delta);
    fprintf(lfp,"max_ix_delta = %ld\n",max_ix_delta);
    fflush(lfp);
  }
  state->ix_delta_sum = ix_delta_sum;
  state->iy_delta_sum = iy_delta_sum;
  state->num_in_tca   = num_in_tca;
}
