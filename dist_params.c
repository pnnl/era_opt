#include "system_includes.h"
#include "state_struct.h"
#include "mpi.h"
#include "unload_params.h"
#include "era_exit.h"
#include "dist_params.h"
int dist_params(struct state_struct *state,
		int filename_length,
		int number_file_names,
		int nproc,
		int  *iparams,
		double *rparams) {
  /*
    Distribute the filenames, the integer and real parameters,
    and unload them into state.
    Called by: era_opt
    Calls:     MPI_Bcast, era_exit, unload_params, fprintf,fflush

    Returns 1 on succes, 0 on failure. 

    Uses the following fields of state:
      output_filename (on rank 0)

    Sets the following fields of state:
      output_filename (on ranks > 0) 
      num_procs           
      num_threads 	     
      source_rows         
      source_columns	     
      channels_rows       
      channels_columns    
      first_channel_index 
      last_channel_index  
      num_channels        
      num_huebins         
      buffer_length       
      num_samples         
      wavelength_index    
      xbar_index	     
      ybar_index          
      zbar_index	     
      xbar_10_index	     
      ybar_10_index	     
      zbar_10_index	     
      s_0_index 	     
      s_1_index	     
      s_2_index	     
      num_radiator_temps  
      max_num_channels    
      nc_stesp_set	     
      num_cct_rows        
      u_target            
      v_target            
      test_rad            
      rf_thresh           
      rad_skosh                 
      watchme
      dump_combos
  */
  int iroot;
  int ierr;

  int success;
  int err_code;

  int rparams_c;
  int iparams_c;

  int gsuccess;
  int fnl;

  int cfnl;
  int agg_fnl;

  fnl       = filename_length;
  success   = 1;
  agg_fnl   = filename_length * number_file_names;
  iroot     = 0;
  /* 
     This works because in alloc0 all of the filename and directory fields are
     allocated consecutively with parameter_file being the first.
     All ranks will have previously called alloc0 to allocate these fields.
  */
  ierr      = MPI_Bcast(state->parameter_file,agg_fnl,MPI_BYTE,iroot,MPI_COMM_WORLD);
  /*
    broadcast all of the file and directory names.
    This works because they are allocated contiguously in alloc0 on all procs.
    Currently the number of these entities is 16. 
    This is hardwired as every proc needs to know it. 
    Might rather have it be an integer parameter that gets distributed first,
    Yet its only going to change when the code base changes, so for now 
    we leave it hardwired.
  */


  /*
    broadcast the output filename base.
  */
  /*
  ierr      = MPI_Bcast(state->output_file,fnl,MPI_BYTE,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    err_code = 50;
    // lfp has not been opened yet. 
    fprintf(stderr,"dist_params, MPI_Bcast of output filenamefailed\n");
    fflush(stderr);
    era_exit(nproc,err_code);
  }
  */
  /*
    broadcast the channels filename base.
  */
  /*
  ierr      = MPI_Bcast(state->channels_file,fnl,MPI_BYTE,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    err_code = 53;
    fprintf(stderr,"dist_params, MPI_Bcast of channels filenamefailed\n");
    fflush(stderr);
    era_exit(nproc,err_code);
  }
  */
  /*
    broadcast the log file filename base.
  */
  /*
  ierr      = MPI_Bcast(state->log_file,fnl,MPI_BYTE,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    err_code = 54;
    fprintf(stderr,"dist_params, MPI_Bcast of log file filename failed\n");
    fflush(stderr);
    era_exit(nproc,err_code);
  }
  */
  /* 
     broadcast the z_dom_file, y_dom_file, and x_dom_file filenames.
  */
  /*
  cfnl = 3 * fnl;
  ierr     = MPI_Bcast(state->z_dom_file,cfnl,MPI_BYTE,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    err_code = 55;
    fprintf(stderr,"dist_params, MPI_Bcast of coordinate file filenames failed\n");
    fflush(stderr);
    era_exit(nproc,err_code);
  }
  */
  rparams_c = 8;
  ierr      = MPI_Bcast(rparams,rparams_c,MPI_DOUBLE,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success = 0;
    err_code = 51;
    fprintf(stderr,"dist_params, MPI_Bcast of rparams failed\n");
    fflush(stderr);
    era_exit(nproc,err_code);
  }
  iparams_c = 64;
  ierr = MPI_Bcast(iparams,iparams_c,MPI_INTEGER,iroot,MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    success  = 0;
    err_code = 52;
    fprintf(stderr,"dist_params, MPI_Bcast of iparams failed\n");
    fflush(stderr);
    era_exit(nproc,err_code);
  } 
  if (success) {
    unload_params(&gsuccess,state,iparams,rparams);
    success = success & gsuccess;
  }
  return (success);
}
