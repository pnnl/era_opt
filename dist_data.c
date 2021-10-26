#include "system_includes.h"
#include "state_struct.h"
#include "mpi.h"
#include "era_exit.h"
#include "dist_data.h"
int dist_data(struct state_struct *state) {
   /*
     Distribute the source_matrix, channels_matrix, huebins, s026_matrix, and
     cct_matrix, isteps_per_channel, subset_indices, and the test_radii_sq  fields of state.
     
     Called by: era_opt
     Calls:     MPI_Bcast,era_exit,fprintf,fflush

     Uses the following fields of state
       channels_rows,
       channels_columns,
       subset_size,
       source_rows,
       source_columns,
       s026_columns,
       num_huebins,
       num_cct_rows,
       use_bruteforce,
       number_scales,

       source_matrix on rank0,
       channel_matrix on rank0,
       subset_indices
       huebins on rank0,
       cct_matrix on rank0,
       s026_matrix on rank0,
       test_radii_sq on rank0,

     Sets the following fields of state on ranks > 0
       *channel_matrix,
       *isteps_per_channel,
       *subset_indices,
       *test_radii_sq,
       *source_matrix,

       *huebins,
       *cct_matrix,
       *s026_matrix,
  */
  int success;
  int iroot;

  int source_rows;
  int source_columns;

  int channels_rows;
  int channels_columns;

  int num_huebins;
  int num_cct_rows;
  
  int bcast_len;
  int ierr;

  int err_code;
  int nproc;

  int num_channels;
  int number_scales;

  int use_multiscale;
  int use_bruteforce;

  int s026_columns;
  int subset_size;
  
  int select_only;
  int ipad;

  FILE *lfp;
  FILE *efp;
  
  success           = 1;
  iroot             = 0;
  nproc             = state->num_procs;
  source_rows       = state->source_rows;
  source_columns    = state->source_columns;
  channels_rows     = state->channels_rows;
  channels_columns  = state->channels_columns;
  s026_columns      = state->s026_columns;
  num_huebins       = state->num_huebins;
  num_cct_rows      = state->num_cct_rows;
  num_channels      = state->num_channels;
  use_multiscale    = state->use_multiscale;
  use_bruteforce    = state->use_bruteforce;
  number_scales     = state->number_scales;
  subset_size       = state->subset_size;
  select_only       = state->select_only;
  lfp               = state->lfp;

  bcast_len        = channels_rows * channels_columns;
  ierr = MPI_Bcast(state->channel_matrix,bcast_len,
		   MPI_DOUBLE,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"dist_data: MPI_Bcast of channel_matrix failed\n");
      fflush(lfp);
    }
    err_code = 72;
    era_exit(nproc,err_code);
  }
  bcast_len = num_channels;
  ierr = MPI_Bcast(state->isteps_per_channel,bcast_len,
		   MPI_INTEGER,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"dist_data_opt: MPI_Bcast of isteps_per_channel failed\n");
      fflush(lfp);
    }
    err_code = 75;
    era_exit(nproc,err_code);
  }
  bcast_len = subset_size;
  ierr = MPI_Bcast(state->subset_indices,bcast_len,
		   MPI_INTEGER,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"dist_data_opt: MPI_Bcast of subset_indices failed\n");
      fflush(lfp);
    }
    err_code = 75;
    era_exit(nproc,err_code);
  }
  bcast_len = num_channels+1;
  ierr = MPI_Bcast(state->test_radii_sq,bcast_len,
		   MPI_DOUBLE,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"dist_data_opt: MPI_Bcast of test_radii_sq failed\n");
      fflush(lfp);
    }
    err_code = 78;
    era_exit(nproc,err_code);
  }
  bcast_len      = source_rows * source_columns;
  ierr = MPI_Bcast(state->source_matrix,bcast_len,
		   MPI_DOUBLE,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"dist_data: MPI_Bcast of source_matrix failed\n");
      fflush(lfp);
    }
    err_code = 71;
    era_exit(nproc,err_code);
  }
  if (select_only == 0) {

    bcast_len        = source_rows * s026_columns;
    ierr = MPI_Bcast(state->s026_matrix,bcast_len,
		     MPI_DOUBLE,iroot,MPI_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"dist_data: MPI_Bcast of s026_matrix failed\n");
	fflush(lfp);
      }
      err_code = 78;
      era_exit(nproc,err_code);
    }
    bcast_len   = 6 * num_huebins;
    ierr = MPI_Bcast(state->huebins,bcast_len,
		     MPI_DOUBLE,iroot,MPI_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"dist_data_opt: MPI_Bcast of huebins failed\n");
	fflush(lfp);
      }
      err_code = 73;
      era_exit(nproc,err_code);
    }
    bcast_len = 3 * num_cct_rows;
    ierr = MPI_Bcast(state->cct_matrix,bcast_len,
		     MPI_DOUBLE,iroot,MPI_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"dist_data_opt: MPI_Bcast of cct_matrix failed\n");
	fflush(lfp);
      }
      err_code = 74;
      era_exit(nproc,err_code);
    }
    /*
    if (use_multiscale == 1) {
      bcast_len = number_scales;
      ierr = MPI_Bcast(state->scales,bcast_len,MPI_DOUBLE,iroot,MPI_COMM_WORLD);
      if (ierr != MPI_SUCCESS) {
        success = 0;
        if (lfp) {
    	fprintf(lfp,"dist_data_opt: MPI_Bcast of scales failed\n");
    	fflush(lfp);
        }
        err_code = 77;
        era_exit(nproc,err_code);
      }
    }
    */
  } /* end if (select_only == 0) */
  return(success);
}
