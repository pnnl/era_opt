#include "system_includes.h"
#include "state_struct.h"
#include "init_z_dom_combos.h"
#include "gen_z_combo.h"
#include "count_nzc.h"
#include "combo_sort.h"
#include "pack_isteps.h"
#include "gen_z_dom_combos.h"
int gen_z_dom_combos(struct state_struct *state) {
  /*
    Called by: gen_combos
    Calls:  init_z_dom_combos, gen_z_combo, count_nzc, 
            combo_sort, pack_isteps
  */
  struct loop_state_struct *loop_state;
  void    *mmap_addr;
  double  *ddata;
  double  *tsv;
  double  *combos;
  double  *sort_scratch;
  double  zcombo[3];
  int64_t *j_ptr;
  int     *istep_index;
  int     *idata;
  int     *imeta_data;
  int64_t mmap_size;
  int64_t j;
  int64_t izc;
  int64_t total_bytes;

  int md_len;
  int bpe;

  int z_combos;
  int mask;
  
  int sbpi;
  int tsv_offset;

  int rec_size;
  int zfd;

  int success;
  int nz;

  int field;
  int flags;

  int ierr;
  int nzc;

  int nz_cpcg_max;
  int ipad;

  FILE *lfp;
  FILE *zfp;
  
  lfp         = state->lfp;
  loop_state  = state->loop_state;
  nz_cpcg_max = state->nz_cpcg_max;
  istep_index = loop_state->istep_index;
  nz          = loop_state->nz;
  /*
    We will use the stored filename in state->z_dom_file
    Then we will mmap this file, its length should be
    ((max_index + 1)^num_z_dom_channels + 4) * record length
    where record length is 20 for now but see thoughts,
    But we add an additional (max_indx+1)^num_z_dom_channels * record_length
    for sorting scratch space.
  */
  success = init_z_dom_combos(state,&ddata,&imeta_data, &zfd, &tsv_offset);
  combos = ddata;
  if (success) {
    /*
      Extract mmap_addr and mmap_size from values set by init_z_dom_combos.
    */
    mmap_addr = (void*)imeta_data;
    j_ptr     = (int64_t*)&imeta_data[6];
    mmap_size = *j_ptr;

    mask = 255;
    sbpi = imeta_data[0] & mask;
    rec_size = tsv_offset + 3;

    izc = (int64_t)0;
    z_combos = 1;
    while (z_combos) {
      z_combos = 0;
      /*
	Loop through z dominant (blue) channels in reverse order (smallest to largest
	z value) setting istep_index
      */
      z_combos = gen_z_combo(loop_state,zcombo);
      /*
	At this point if z_combos = 1
	then loop_state->istep_index[0:nz-1] contains a valid set of steps
	for a combination of z_dom channels, and zcombo has
	the combined xyz contributions of those channels,
	And loop_state->prev_xyz and loop_state->istep_index indicate 
	starting values for the next call to gen_z_combo.c
      */
      if (z_combos) {
	/* 
	   Check for number of nonzero channels per coordinate group
	*/
	nzc = count_nzc(nz,istep_index);
	if (nzc <= nz_cpcg_max) {
	  izc = izc + 1;
	  /*
	    Record istep_index and tristate values for this 
	    combination.
	  */
	  idata = (int*)ddata;
	  tsv = (double *)&ddata[tsv_offset];
	  pack_isteps(nz,sbpi,istep_index,idata);
	  tsv[0] = zcombo[0];
	  tsv[1] = zcombo[1];
	  tsv[2] = zcombo[2];
	  ddata += rec_size; /* Caution address arithmetic */
	}
      }
    }
    sort_scratch = ddata;
    field = tsv_offset + 2;
    combo_sort(izc,rec_size,field,combos,sort_scratch);
    j_ptr = (int64_t*)&imeta_data[8];
    *j_ptr = izc;
    /* 
      Now we should set the length of the color file to its true length. 
      which is md_len + izc * bpe 
    */
    md_len = (imeta_data[0] >> 24) & mask;
    bpe    = (imeta_data[1] >> 8) & mask;
    total_bytes = md_len + (izc * bpe);
    j_ptr = (int64_t*)&imeta_data[6];
    *j_ptr = total_bytes;
    /*
      Push the file out to disk.
    */
    flags = MS_SYNC;
    msync(mmap_addr,mmap_size,flags);
    munmap(mmap_addr,mmap_size);
    fclose(state->zfp);
    /* 
      Truncate the file to its true length (removing the sort scratch).
    */
    ierr = truncate(state->z_dom_file,total_bytes);
    if (ierr != 0) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"gen_z_dom_combos: Error truncationg %s to %ld bytes\n",
		state->z_dom_file,total_bytes);
	fflush(lfp);
      }
    }
    if(lfp) {
      fprintf(lfp,"izc = %ld\n",izc);
      fflush(lfp);
    }
  }
  return(success);
}
