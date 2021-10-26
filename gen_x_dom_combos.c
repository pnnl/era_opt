#include "system_includes.h"
#include "state_struct.h"
#include "init_x_dom_combos.h"
#include "gen_x_combo.h"
#include "count_nzc.h"
#include "combo_sort.h"
#include "pack_isteps.h"
#include "gen_x_dom_combos.h"
int gen_x_dom_combos(struct state_struct *state) {
  /*
    Generate all possible combinations of x-dominant tsv values
    at all power levels.
    Called by: gen_combos
    Calls:     init_x_dom_combos, gen_c_combo, count_nzc, 
               combo_sort, pack_isteps
  */
  struct loop_state_struct *loop_state;
  void    *mmap_addr;
  double  *ddata;
  double  *tsv;
  double  *combos;
  double  *sort_scratch;
  double  *xvals;
  double  *r_ptr;
  double  xcombo[5];
  int64_t *xbins;
  int64_t *j_ptr;
  int     *istep_index;
  int     *idata;
  int     *imeta_data;
  double  xtop;
  double  xbot;
  double  x_dom_vb_min;
  double  x_dom_vb_max;
  double  average_combos_per_bin;
  double  sd;
  double  dev;
  double  combo_bin_width;
  double  bin_top;
  double  xmax;
  int64_t mmap_size;
  int64_t j;
  int64_t ixc;
  int64_t bin_ct;
  int64_t total_bytes;

  int x_combos;
  int mask;
  
  int sbpi;
  int tsv_offset;

  int nfields;
  int xfd;

  int success;
  int nz;

  int field;
  int flags;

  int num_bins;
  int bin_pos;

  int i;
  int max_combos_per_bin;

  int combos_in_bin;
  int nx;

  int ny;
  int md_len;

  int ierr;
  int bpe;

  int nz_cpcg_max;
  int nzc;


  FILE *lfp;
  FILE *zfp;
  
  lfp         = state->lfp;
  nz_cpcg_max = state->nz_cpcg_max;
  loop_state  = state->loop_state;
  istep_index = loop_state->istep_index;
  nz          = loop_state->nz;
  ny          = loop_state->ny;
  nx          = loop_state->nx;
  xmax        = loop_state->x_top;
  /*
    We need to create a file name for the blue channel combos
    Extract the channel filename base. add a .z.dat extension
    for the output filename.
    Then we will mmap this file, its length should be
    ((max_index + 1)^num_blue_channels + 4) * record length
    where record length is 20 for now but see thoughts,
    But we add an additional (max_indx+1)^num_blue_channels * record_length
    for sorting scratch space.
  */
  success = init_x_dom_combos(state,&ddata,&imeta_data, &xfd, &tsv_offset);
  combos = ddata;
  if (success) {
    /*
      Extract mmap_addr and mmap_size from values set by init_blue_combos.
    */
    mmap_addr = (void*)imeta_data;
    j_ptr     = (int64_t*)&imeta_data[6];
    mmap_size = *j_ptr;
    num_bins = imeta_data[5];

    mask = 255;
    sbpi = imeta_data[0] & mask;
    nfields = tsv_offset + 5;

    ixc = (int64_t)0;
    x_combos = 1;
    while (x_combos) {
      x_combos = 0;
      /*
	Loop through z (blue) channels in reverse order (smallest to largest
	z value) setting istep_index
      */
      x_combos = gen_x_combo(loop_state,xcombo);
      /*
	At this point if y_combos = 1
	then loop_state->istep_index[nz:nz+ny-1] contains a valid set of steps
	for a combination of green channels, and xcombo has
	the combined xyz contributions of those channels,
	And loop_state->prev_xyz and loop_state->istep_index indicate 
	starting values for the next call to gen_y_combo.c
      */
      if (x_combos) {
	/* 
	   Check for number of nonzero channels per coordinate group
	*/
	nzc = count_nzc(nx,&istep_index[ny+nz]);
	if (nzc <= nz_cpcg_max) {
	  ixc = ixc + 1;
	  /*
	    Record istep_index and tristate values for this 
	    combination.
	  */
	  idata = (int*)ddata;
	  tsv = (double *)&ddata[tsv_offset];
	  pack_isteps(nx,sbpi,&istep_index[ny+nz],idata);
	  tsv[0] = xcombo[0];
	  tsv[1] = xcombo[1];
	  tsv[2] = xcombo[2];
	  tsv[3] = xcombo[3];
	  tsv[4] = xcombo[4];
	  ddata += nfields; /* Caution address arithmetic */
	}
      }
    } /* end while y combos */
    sort_scratch = ddata;
    field = tsv_offset+4;
    combo_sort(ixc,nfields,field,combos,sort_scratch);
    j_ptr = (int64_t*)&imeta_data[8];
    *j_ptr = ixc;
    xtop = combos[(ixc-1)*nfields + field];
    xbot = combos[field];
    /*
      make xbot an even multiple of combo_bin_width;
    xbot = ((double)((int)(xbot/combo_bin_width)))*combo_bin_width;
    num_bins = (int)(((xtop-xbot)/combo_bin_width) + 1.5);
    num_bins now computed in init_x_dom_combos.
    */
    combo_bin_width = (xtop - xbot)/(num_bins-1);
    /*
      We need to store the combo_bin_width in the meta data.
    */
    r_ptr = (double*)&imeta_data[14];
    *r_ptr = combo_bin_width;
    bin_pos = 0;
    bin_ct  = 0;
    xbins   = (int64_t*)sort_scratch;
    xvals   = &combos[field];
    bin_top = combo_bin_width;
    xbins[0] = (int64_t)0;
    for (i=1;i<num_bins;i++) {
      bin_top = (i * combo_bin_width) + xbot;
      while ((*xvals < bin_top) && (bin_ct < ixc)) {
	xvals += nfields;
	bin_ct += (int64_t)1;
      }
      xbins[i] = bin_ct;
    }
    xbins[num_bins - 1] = ixc;
    /*
      Store the offeset in 8byte wores of ybins from combos (ddata)
      in the meta data.
    */
    j_ptr = (int64_t*)&imeta_data[10];
    *j_ptr = nfields * ixc;
    /*
      Reset the total bytes in the meta data field
      total bytes = md_len + ixc * bpe + num_bins * 8
      md_len = (imeta_data[0] >> 24) & 255
      bpe    = (imeta_data[1] >> 8) & 255  
    */
    mask   = 255;
    md_len = (imeta_data[0] >> 24) & mask;
    bpe    = (imeta_data[1] >> 8) & mask;
    total_bytes = md_len + (ixc * bpe) + (num_bins * 8);
    j_ptr = (int64_t*)&imeta_data[6];
    *j_ptr = total_bytes;
    /*
      Now we want to gather some statistics about y bins.
      Specifically we want to know the maximum number of combos in a bin,
      the average number of combos in a bin and possibly the 
      standard deviation of the number of combos in a bin.
    */
    max_combos_per_bin = 0;
    average_combos_per_bin = ((double)ixc)/((double)num_bins);
    sd  = 0;
    for (i=0;i<num_bins-2;i++) {
      combos_in_bin = xbins[i+1] - xbins[i];
      if (combos_in_bin > max_combos_per_bin) {
	max_combos_per_bin = combos_in_bin;
      }
      dev = (combos_in_bin - average_combos_per_bin);
      sd = sd + (dev*dev);
    }
    sd = sqrt(sd/(num_bins-1));

    /*
      Now we need to store ybot, and ytop as the green_vb_min and
      green_vb_max values in the meta data.
    */
    r_ptr = (double*)&imeta_data[20];
    *r_ptr = xbot;
    r_ptr += 1; /* caution address arithmetic */
    *r_ptr = xtop;
    /*
      We also want the maximum ub_val and minimum ub_val for these
      combinations:
    */
    field        = tsv_offset + 3;
    xvals        = &combos[field];
    x_dom_vb_min = *xvals;
    x_dom_vb_max = *xvals;
    xvals += nfields; /* Caution address arithmetic */
    for (i=1;i<ixc;i++) {
      if (*xvals < x_dom_vb_min) {
	x_dom_vb_min = *xvals;
      } else {
	if (*xvals > x_dom_vb_max) {
	  x_dom_vb_max = *xvals;
	}
      }
      xvals += nfields; /* Caution address arithmetic */
    }
    r_ptr = (double*)&imeta_data[16];
    *r_ptr = x_dom_vb_min;
    r_ptr += 1; /* caution address arithmetic */
    *r_ptr = x_dom_vb_max;

    flags = MS_SYNC;
    msync(mmap_addr,mmap_size,flags);
    munmap(mmap_addr,mmap_size);
    fclose(state->xfp);
    /* 
      Truncate the file to its true length (removing the excess sort scratch).
    */
    ierr = truncate(state->x_dom_file,total_bytes);
    if (ierr != 0) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"gen_x_dom_combos: Error truncationg %s to %ld bytes\n",
		state->x_dom_file,total_bytes);
	fflush(lfp);
      }
    }
    if(lfp) {
      fprintf(lfp,"ixc = %ld\n",ixc);
      fprintf(lfp,"num x_dom bins = %d\n",num_bins);
      fprintf(lfp,"xtop (x_dom_ub_max) = %le, xbot (x_dom_ub_min) = %le\n",
	      xtop,xbot);
      fprintf(lfp,"x_dom_vb_min = %le\n",x_dom_vb_min);
      fprintf(lfp,"x_dom_vb_max = %le\n",x_dom_vb_max);
      fprintf(lfp,"average_combos_per_x_dom_bin = %le\n",
	      average_combos_per_bin);
      fprintf(lfp,"max_combos_per_x_dom_bin     = %d\n",
	      max_combos_per_bin);
      fprintf(lfp,"sd of combos_per_x_dom_bin   = %le\n",sd);
      fflush(lfp);
    }
  }
  return(success);
}
