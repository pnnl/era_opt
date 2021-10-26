#include "system_includes.h"
#include "state_struct.h"
#include "unpack_isteps.h"
#include "count_nzc.h"
#include "process_combo.h"
/*
#include "unique_ratio.h"
#include "get_steps.h"
#include "get_ipsteps.h"
#include "ies_tm_30_18.h"
#include "output_tm_30_results.h"
#include "output_tm_30_debug.h"
#include "step_compare.h"
*/
#include "eprocess_combos.h"
void eprocess_combos(struct state_struct *state,
		     struct loop_state_struct *loop_state, 
		     int myid, int nproc) {
  /*
    Triply nested set of pnnl loops to test all combinations.
    Called by: era_opt
    Calls:     count_nzc, unpack_isteps, unique_ratio, get_ipsteps, 
    Also calls get_steps, ies_tm_30_18, output_tm_30_results,       
               output_tm_30_debug, step_compare when select_only is
	       not set. 
  */
  /*
  double *stepped_channel_xyz;
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
  double *vb_plane;
  double *ub_plane;
  double *dbg_steps;
  double *pstep_fracts;
  double *stest;
  double *channels;
  double cct_vals[4];

  int64_t *iy_dom_bins; /* state field */
  int64_t *ix_dom_bins;   /* state field */

  int *istep_index;
  int *ipstep_index;
  int *mm_istep_index;
  int *imap_sort;
  int *unique_prime_factors;
  int *upf_index;
  int *iz_dom_steps;
  int *iy_dom_steps;
  int *ix_dom_steps;

  double ymax;
  double xmax;
  double ymin;
  double xmin;
  double ymin_c;
  double xmin_c;
  double ymax_c;
  double xmax_c;
  double z;
  double y;
  double x;
  double uprime;
  double vprime;
  double denom;
  double rdenom;
  double z_dom_vb;
  double z_dom_ub;
  double y_dom_ub;
  double teal_ub;
  double y_dom_combo_bin_width;
  double y_dom_vb_min;
  double x_dom_combo_bin_width;
  double x_dom_vb_max;
  double x_dom_vb_min;
  double x_dom_ub_min;
  double ua_delta;
  double va_delta;
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
  double thresh;

  int64_t izc;
  /*
  int64_t iyc;
  int64_t ixc;
  */
  int64_t iy_start;
  int64_t iy_end;
  int64_t ix_start;
  int64_t ix_end;

  int64_t iz;
  int64_t iy;
  int64_t ix;

  int64_t iy_delta_sum;
  int64_t ix_delta_sum;
  int64_t iy_delta;
  int64_t ix_delta;
  int64_t num_in_tca;
  int64_t lix_delta_sum;
  int64_t max_ix_delta;
  int64_t max_iy_delta;

  int num_channels;
  int nz;

  int ny;
  /*
  int z_combos;

  int y_combos;
  int x_combos;
  */
  /*
  int istep;
  int i3step;
  */
  /*
  int ithree;
  int inc1;
  */
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
  
  /*
  int iymin_bin;
  int iymax_bin;
  int ixmin_bin;
  int ixmax_bin;
  */

  int nzzc;
  int nzyc;
  
  int nztc;
  int nzxc;

  int unique;
  int max_index;

  int full_output;
  int source_rows;

  int j;
  int dump_combos;

  /*
  char n_char;
  */
  FILE *lfp;
  FILE *ofp;
  FILE *dfp;
  FILE *efp;

  source_rows         = state->source_rows;
  istep_index  	      = loop_state->istep_index;
  ipstep_index        = loop_state->ipstep_index;
  mm_istep_index      = loop_state->mm_istep_index;
  imap_sort           = loop_state->imap_sort;
  unique_prime_factors = loop_state->unique_prime_factors;
  upf_index           = loop_state->upf_index;
  max_index           = loop_state->max_index;
  num_channels 	      = loop_state->num_channels;
  nz                  = loop_state->nz;
  ny                  = loop_state->ny;
  nx                  = loop_state->nx;
  xmin_c              = loop_state->xmin_c;
  ymin_c              = loop_state->ymin_c;
  xmax_c              = loop_state->xmax_c;
  ymax_c              = loop_state->ymax_c;
  lfp                 = loop_state->lfp;
  fracts              = loop_state->fracts;
  va_delta            = loop_state->va_delta;
  ua_delta            = loop_state->ua_delta;
  pstep_fracts        = loop_state->pstep_fracts;

  tm_30_results       = state->tm_30_results;
  rf_thresh           = state->rf_thresh;
  vb_plane            = &loop_state->vb_plane[0];
  ub_plane            = &loop_state->ub_plane[0];
  ofp                 = state->ofp;
  extra_log_info      = state->extra_log_info;
  test_radii_sq       = state->test_radii_sq;
  ybar_t_channels     = state->ybar_t_channels;
  e_t_channels        = state->e_t_channels;
  sz_t_channels       = state->sz_t_channels;
  full_output         = state->full_output;
  dbg_steps           = state->dbg_steps;
  izc = 0;
  u_target = loop_state->u_target;
  v_target = loop_state->v_target;
  
  iz_dom_sbpi        = state->iz_dom_sbpi;
  iz_dom_tsv_offset  = state->iz_dom_tsv_offset;
  nz_dom_fields      = state->nz_dom_fields;
  izc                = state->izc;
  dump_combos        = state->dump_combos;
  dfp                = state->dfp;
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
  y_dom_combo_bin_width = state->y_dom_combo_bin_width;
  y_dom_vb_min      = state->y_dom_vb_min;
                   
  ix_dom_sbpi         = state->ix_dom_sbpi;
  ix_dom_tsv_offset   = state->ix_dom_tsv_offset;
  nx_dom_fields       = state->nx_dom_fields;
  nx_dom_bins         = state->nx_dom_bins;

  x_dom_data          = state->x_dom_data;
  ix_dom_bins         = state->ix_dom_bins;
  x_dom_combo_bin_width = state->x_dom_combo_bin_width;
  x_dom_vb_min        = state->x_dom_vb_min;
  x_dom_vb_max        = state->x_dom_vb_max;
  x_dom_ub_min        = state->x_dom_ub_min;
  watchme             = state->watchme;

  num_in_tca   = (int64_t)0;
  iy_delta_sum = (int64_t)0;
  ix_delta_sum = (int64_t)0;
  max_iy_delta = (int64_t)0;
  max_ix_delta = (int64_t)0;
  iz_dom_delta = nz_dom_fields * nproc;
  my_z_dom_data = z_dom_data + (myid * nz_dom_fields); /* Caution, Address arithmetic */
  progress = 0.0;
  for (iz = myid; iz < izc; iz += nproc) {
    iz_dom_steps = (int*)my_z_dom_data;
    unpack_isteps(nz,iz_dom_sbpi,istep_index,iz_dom_steps);
    nzzc   = count_nzc(nz,istep_index);
    z_dom_x = my_z_dom_data[iz_dom_tsv_offset];
    z_dom_y = my_z_dom_data[iz_dom_tsv_offset+1];
    z_dom_z = my_z_dom_data[iz_dom_tsv_offset+2];
    z_dom_vb = (z_dom_x * vb_plane[0]) + 
              (z_dom_y * vb_plane[1]) + 
              (z_dom_z * vb_plane[2]);
    z_dom_ub = (z_dom_x * ub_plane[0]) + 
              (z_dom_y * ub_plane[1]) + 
              (z_dom_z * ub_plane[2]);
    /*
      Round robin scheduling of tasks.
    */
    my_z_dom_data += iz_dom_delta; /* Caution, Address arithmetic. */
    /*
      Determine  the y_dom_vb range. We want vb value to be < 0 but > -dist((xmax,ymax),(xc,yc)) at zmax

      We also have to remember that the red channel can contribute between x_dom_vb_min and x_dom_vb_max to the
      vb value.
      So we want         vb = z_dom_vb + y_dom_vb + x_dom_vb to be in [-dist((xmax,ymax),(xc,yc)) at zmax,0]


                                  we have 

	  z_dom_vb + y_dom_vb + x_dom_min_vb < z_dom_vb+y_dom_vb+x_dom_vb < z_dom_vb+y_dom_vb+x_dom_vb_max
         
           with x_dom_vb in [x_dom_vb_min,x_dom_vb_max]
           let va_delta = -dist((xmax,ymax),(xc,yc)) at zmax
         
	   so we want y_dom_vb in [-va_delta - z_dom_vb - x_dom_vb_max < y_dom_vb < -z_dom_vb - x_dom_vb_min]

	   now we need to understand which bin to start in and which to end in.

	   Now the green channels are binned from y_dom_vb_min to y_dom_vb_max in increments of combo_bin_width

	   ivb_low = (int)((-va_delta - z_dom_vb - x_dom_vb_max - y_dom_vb_min)/combo_bin_width)

	   and ivb_high is (int)(((-z_dom_vb - x_dom_vb_min- y_dom_vb_min)/combo_bin_width) + 1.)
	   
	   and we need vb_low >= 0 and vb_high < ny_dom_bins;
         
    */
    ivb_low = (int)((-va_delta - z_dom_vb - x_dom_vb_max - y_dom_vb_min)/y_dom_combo_bin_width);
    ivb_high = (int)(((-z_dom_vb -x_dom_vb_min - y_dom_vb_min)/y_dom_combo_bin_width) + 1.0);
    if (ivb_low < 0) {
      ivb_low = 0;
    }
    if (ivb_low > ny_dom_bins-1) {
      ivb_low = ny_dom_bins-1;
    }
    if (ivb_high > (ny_dom_bins-1)) {
      ivb_high = ny_dom_bins-1;
    }
    if (ivb_high < 0) {
      ivb_high = 0;
    }
    iy_start = iy_dom_bins[ivb_low];
    iy_end   = iy_dom_bins[ivb_high];
    iy_delta = iy_end - iy_start;
    if (iy_delta > max_iy_delta) {
      max_iy_delta = iy_delta;
    }
    iy_delta_sum += iy_delta;
    ymax = ymax_c * z_dom_z;
    ymin = ymin_c * z_dom_z;
    /*
      Subtract off the existing y.
    */
    ymax = ymax - z_dom_y;
    ymin = ymin - z_dom_y;
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
      /*
      y_dom_vb = my_y_dom_data[iy_dom_tsv_offset+3];
      */
      y_dom_ub = my_y_dom_data[iy_dom_tsv_offset+4];
      teal_x = z_dom_x + y_dom_x;
      teal_y = z_dom_y + y_dom_y;
      teal_z = z_dom_z + y_dom_z;
      teal_ub = z_dom_ub + y_dom_ub;
      my_y_dom_data += ny_dom_fields; /* Caution Address arithmetic */
      /*
	Now for this blue and green combination we want to examine 
	all of the red combinataions that land above line ub and below
	line ub+ua_delta;

	Thus we want teal_ub + x_dom_ub to be in [0,ua_delta]
	Hence we want red combos with x_dom_ub in [-teal_ub,ua_delta-teal_ub]

	iub_low = (int)((-teal_ub-x_dom_ub_min)/combo_bin_width);

      */
      iub_low = (int)((-teal_ub-x_dom_ub_min)/x_dom_combo_bin_width);
      iub_high = (int)(((ua_delta-teal_ub-x_dom_ub_min)/x_dom_combo_bin_width) + 1.0);
      if (iub_low < 0) {
	iub_low = 0;
      }
      if (iub_low > nx_dom_bins-1) {
	iub_low = nx_dom_bins-1;
      }
      if (iub_high > (nx_dom_bins-1)) {
	iub_high = nx_dom_bins - 1;
      }
      if (iub_high < 0) {
	iub_high = 0;
      }
      /*
#define DBG 1
      */
      /*
#ifdef DBG
      fprintf(lfp,"iub_low = %d,iub_high=%d\n",iub_low,iub_high);
      fflush(lfp);
#endif
      */
      xmax = xmax_c * teal_z;
      xmin = xmin_c * teal_z;
      xmax = xmax - teal_x;
      xmin = xmin - teal_x;
      /*
      ixmin_bin = (xmin/combo_bin_width);
      ixmax_bin = ((xmax/combo_bin_width) + .5);
      ix_start = ix_dom_bins[ixmin_bin];
      ix_end   = ix_dom_bins[ixmax_bin];
      if (ix_end >= ixc) {
	ix_end = ixc - 1;
      }
      */
      ix_start = ix_dom_bins[iub_low];
      ix_end   = ix_dom_bins[iub_high];
      ix_delta = ix_end - ix_start + 1;
      lix_delta_sum += ix_delta;
      my_x_dom_data = x_dom_data + (ix_start * nx_dom_fields); /* Cautions Addersss arithmetic */
      /*
#ifdef DBG
      fprintf(lfp,"ix_start = %ld, ix_end = %ld\n",ix_start,ix_end);
      fflush(lfp);
#endif
      */
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
	/*
	  Now we have an x,y,z candiate with the power steps that generated
	  it in istep_index.
	*/
	/* 
	   Not sure we need to do the following
	   But it possibly recomputes x,y,z more accurately for these steps.
	   for (i=0;i<num_channels;i++) {
  	     step_fracts[i] = fracts[i] * istep_index[i];
           }
	   dgemv_(&n_char,&ithree,&num_channels,
	          &alpha,xyz_sort,&ithree,
                  step_fracts,&inc1,
    	          &beta,curr_xyz,&inc1,1);
           x = curr_xyz[0];
	   y = curr_xyz[1];
	   z = curr_xyz[2];
	   End of x,y,z recalc.
	*/
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
	  /*
	  unique = unique_ratio(num_channels,max_index,istep_index,
				unique_prime_factors,upf_index,mm_istep_index);

	  if (unique) {
	    num_in_tca += (int64_t)1;
	    steps = loop_state->step_fracts;
	    get_steps(num_channels,imap_sort,mm_istep_index,fracts,steps);
#ifdef DBG
	    state->print_efit = 0;
	    thresh = 0.001;
	    if (step_compare(num_channels,thresh,dbg_steps,steps) == 1) {
	      fprintf(lfp,"Debug check");
	      state->print_efit = 1;
	    }
	    if (state->print_efit) {
	      stest = state->stest;
	      for (i=0;i<source_rows;i++) {
		stest[i] = 0.0;
	      }
	      channels = state->channels;
	      for (j = 0;j<num_channels;j++) {
		for (i=0;i<source_rows;i++) {
		  stest[i] += steps[j]*channels[i];
		}
		channels = channels + source_rows; // Caution address arithmetic 
	      }
	      fprintf(lfp,"stest = \n");
	      for (i=0;i<source_rows;i++) {
		fprintf(lfp,"stest[%d] = %le\n",
			i,stest[i]);
	      }
	      fflush(lfp);
	    }
#endif
	    ies_tm_30_18(state,steps,tm_30_results,cct_vals);
	    if (tm_30_results[0] >= rf_thresh) {
	      if (full_output) {
		ybar_dot_stest = ybar_t_channels[0] * steps[0];
		for(i=1;i<num_channels;i++) {
		  ybar_dot_stest += ybar_t_channels[i] * steps[i];
		}
		stest_sum      = e_t_channels[0] * steps[0];
		for (i=1;i<num_channels;i++) {
		  stest_sum    += e_t_channels[i] * steps[i];
		}
		ler   = 683.0*ybar_dot_stest/stest_sum;
		sum_m = sz_t_channels[0] * steps[0];
		for (i=1;i<num_channels;i++) {
		  sum_m += sz_t_channels[i] * steps[i];
		}
		sum_p = ybar_dot_stest;
		m_over_p = sum_m/sum_p;
		dist = sqrt(sqdist);
		output_tm_30_results(state,loop_state,
				     steps,tm_30_results,ler,m_over_p,sum_p,dist,
				     cct_vals);
		if (extra_log_info) {
		  get_ipsteps(num_channels,imap_sort,mm_istep_index,ipstep_index);
		  output_tm_30_debug(state,lfp,num_channels,ipstep_index,
				     dist,tm_30_results,nzc,ler,m_over_p,sum_p);
		}
	      }
	    } // end if rf_thresh is passed 
	  } // end if unique 
	  */
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
