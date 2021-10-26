#include "system_includes.h"
#include "state_struct.h"
#include "init_y_dom_combos.h"
#include "gen_y_combo.h"
#include "count_nzc.h"
#include "combo_sort.h"
#include "pack_isteps.h"
#include "gen_y_dom_combos.h"
int gen_y_dom_combos(struct state_struct *state) {
  struct loop_state_struct *loop_state;
  void    *mmap_addr;
  double  *ddata;
  double  *tsv;
  double  *combos;
  double  *sort_scratch;
  double  *yvals;
  double  *r_ptr;
  double  ycombo[5];
  int64_t *ybins;
  int64_t *j_ptr;
  int     *istep_index;
  int     *idata;
  int     *imeta_data;
  double  ytop;
  double  ybot;
  double  ymax;
  double  average_combos_per_bin;
  double  sd;
  double  dev;
  double  combo_bin_width;
  double  bin_top;
  double  y_dom_ub_min;
  double  y_dom_ub_max;
  int64_t mmap_size;
  int64_t j;
  int64_t iyc;
  int64_t bin_ct;
  int64_t total_bytes;

  int y_combos;
  int mask;
  
  int sbpi;
  int tsv_offset;

  int nfields;
  int yfd;

  int nz;
  int ny;

  int field;
  int flags;

  int num_bins;
  int bin_pos;

  int i;
  int max_combos_per_bin;

  int combos_in_bin;
  int md_len;

  int bpe;
  int ierr;

  int success;
  int nz_cpcg_max;
  
  int nzc;
  int ipad;


  FILE *lfp;
  FILE *yfp;
  
  lfp         = state->lfp;
  loop_state  = state->loop_state;
  nz_cpcg_max = state->nz_cpcg_max;
  istep_index = loop_state->istep_index;
  nz          = loop_state->nz;
  ny          = loop_state->ny;
  ymax        = loop_state->y_top;
  /*
    We will use the file name in state->y_dom_file for a filename.
    Then we will mmap this file, its length should be
    ((max_index + 1)^num_blue_channels + 4) * record length
    where record length is 20 for now but see thoughts,
    But we add an additional (max_indx+1)^num_y_dom_channels * record_length
    for sorting scratch space, but we also need to know that we have 
    enough space for the bins storage which is num_bins * 8 where 
    num_bins =ceil((ymax/combo_bin_width)+1.5)
  */
  success = init_y_dom_combos(state,&ddata,&imeta_data, &yfd, &tsv_offset);
  combos = ddata;
  if (success) {
    /*
      Extract mmap_addr and mmap_size from values set by init_y_dom_combos.
    */
    mmap_addr = (void*)imeta_data;
    j_ptr     = (int64_t*)&imeta_data[6];
    mmap_size = *j_ptr;

    mask = 255;
    sbpi = imeta_data[0] & mask;
    nfields = tsv_offset + 5;
    num_bins = imeta_data[5];
    iyc = (int64_t)0;
    y_combos = 1;
    while (y_combos) {
      y_combos = 0;
      /*
	Loop through y-dominant (green) channels in reverse order (smallest to largest
	z value) setting istep_index
      */
      y_combos = gen_y_combo(loop_state,ycombo);
      /*
	At this point if y_combos = 1
	then loop_state->istep_index[nz:nz+ny-1] contains a valid set of steps
	for a combination of y-dominant channels, and ycombo has
	the combined xyz contributions of those channels,
	And loop_state->prev_xyz and loop_state->istep_index indicate 
	starting values for the next call to gen_y_combo.c
      */
      if (y_combos) {
	/* 
	   Check for number of nonzero channels per coordinate group
	*/
	nzc = count_nzc(ny,&istep_index[nz]);
	if (nzc <= nz_cpcg_max) {
	  
	  iyc = iyc + 1;
	  /*
	    Record istep_index and tristate values for this 
	    combination.
	  */
	  idata = (int*)ddata;
	  tsv = (double *)&ddata[tsv_offset];
	  pack_isteps(ny,sbpi,&istep_index[nz],idata);
	  tsv[0] = ycombo[0];
	  tsv[1] = ycombo[1];
	  tsv[2] = ycombo[2];
	  tsv[3] = ycombo[3];
	  tsv[4] = ycombo[4];
	  ddata += nfields; /* Caution address arithmetic */
	}
      }
    } /* end while y combos */
    /*
      y_combo data is in x,y,z,vb_value,ub_value order
      Sort on vb_value
    */
    sort_scratch = ddata;
    field = tsv_offset + 3;
    combo_sort(iyc,nfields,field,combos,sort_scratch);
    j_ptr = (int64_t*)&imeta_data[8];
    *j_ptr = iyc;
    ytop = combos[(iyc-1)*nfields + field];
    ybot = combos[field];
    combo_bin_width = (ytop - ybot)/(num_bins-1);
    /*
      We need to store the combo_bin_width in the meta data.
    */
    r_ptr = (double*)&imeta_data[14];
    *r_ptr = combo_bin_width;
    bin_pos = 0;
    bin_ct  = 0;
    ybins   = (int64_t*)sort_scratch;
    yvals   = &combos[field];
    ybins[0] = (int64_t)0;
    for (i=1;i<num_bins;i++) {
      bin_top = (i * combo_bin_width) + ybot;
      while ((*yvals < bin_top) && (bin_ct < iyc)) {
	yvals += nfields;
	bin_ct += (int64_t)1;
      }
      ybins[i] = bin_ct;
    }
    ybins[num_bins-1] = iyc;
    /*
      Store the offeset in 8byte words of ybins from combos (ddata)
      in the meta data.
    */
    j_ptr = (int64_t*)&imeta_data[10];
    *j_ptr = nfields * iyc;
    /*
      Reset the total bytes in the meta data field
      total bytes = md_len + iyc * bpe + num_bins * 8
      md_len = (imeta_data[0] >> 24) & 255
      bpe    = (imeta_data[1] >> 8) & 255  
    */
    mask   = 255;
    md_len = (imeta_data[0] >> 24) & mask;
    bpe    = (imeta_data[1] >> 8) & mask;
    total_bytes = md_len + (iyc * bpe) + (num_bins * 8);
    j_ptr = (int64_t*)&imeta_data[6];
    *j_ptr = total_bytes;
    /*
      Now we want to gather some statistics about y bins.
      Specifically we want to know the maximum number of combos in a bin,
      the average number of combos in a bin and possibly the 
      standard deviation of the number of combos in a bin.
    */
    max_combos_per_bin = 0;
    average_combos_per_bin = ((double)iyc)/((double)num_bins);
    sd  = 0.0;
    for (i=0;i<num_bins-2;i++) {
      combos_in_bin = ybins[i+1] - ybins[i];
      if (combos_in_bin > max_combos_per_bin) {
	max_combos_per_bin = combos_in_bin;
      }
      dev = (combos_in_bin - average_combos_per_bin);
      sd = sd + (dev*dev);
    }
    sd = sqrt(sd/(num_bins-1));
    /*
      Now we need to store ybot, and ytop as the y_dom_vb_min and
      y_dom__vb_max values in the meta data.
    */
    r_ptr = (double*)&imeta_data[16];
    *r_ptr = ybot;
    r_ptr += 1; /* caution address arithmetic */
    *r_ptr = ytop;
    /*
      We also want the maximum ub_val and minimum ub_val for these
      combinations:
    */
    field        = tsv_offset + 4;
    yvals        = &combos[field];
    y_dom_ub_min = *yvals;
    y_dom_ub_max = *yvals;
    yvals += nfields; /* Caution address arithmetic */
    for (i=1;i<iyc;i++) {
      if (*yvals < y_dom_ub_min) {
	y_dom_ub_min = *yvals;
      } else {
	if (*yvals > y_dom_ub_max) {
	  y_dom_ub_max = *yvals;
	}
      }
      yvals += nfields; /* Caution address arithmetic */
    }
    r_ptr += 1; /* caution address arithmetic */
    *r_ptr = y_dom_ub_min;
    r_ptr += 1; /* caution address arithmetic */
    *r_ptr = y_dom_ub_max;

    flags = MS_SYNC;
    msync(mmap_addr,mmap_size,flags);
    munmap(mmap_addr,mmap_size);
    fclose(state->yfp);
    /* 
      Truncate the file to its true length (removing the excess sort scratch).
    */
    ierr = truncate(state->y_dom_file,total_bytes);
    if (ierr != 0) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"gen_y_dom_combos: Error truncationg %s to %ld bytes\n",
		state->y_dom_file,total_bytes);
	fflush(lfp);
      }
    }
    if(lfp) {
      fprintf(lfp,"iyc = %ld\n",iyc);
      fprintf(lfp,"num y_dom bins = %d\n",num_bins);
      fprintf(lfp,"ytop (y_dom_vb_max) = %le, ybot (y_dom_vb_min) = %le\n",
	      ytop,ybot);
      fprintf(lfp,"y_dom_ub_min = %le\n",y_dom_ub_min);
      fprintf(lfp,"y_dom_ub_max = %le\n",y_dom_ub_max);
      fprintf(lfp,"average_combos_per_y_dom_bin = %le\n",
	      average_combos_per_bin);
      fprintf(lfp,"max_combos_per_y_dom_bin     = %d\n",
	      max_combos_per_bin);
      fprintf(lfp,"sd of combos_per_y_dom_bin   = %le\n",sd);
      fflush(lfp);
    }
  }
  return(success);
}
