#include "system_includes.h"
#include "state_struct.h"
#include "alloc0.h"
#include "read_params.h"
#include "alloc1.h"
#include "fill_isteps_per_channel.h"
#include "read_data.h"
#include "build_coord_filenames.h"
#include "prep.h"
/*
#include "write_csv.h"
*/
#include  "init_loop_state.h"
/*
#include  "gen_blue_combos.h"
#include  "gen_green_combos.h"
#include  "gen_red_combos.h"
*/
#include "bf_gen_z_dom_combos.h"
#include "gen_y_dom_combos.h"
#include "gen_x_dom_combos.h"
#include "echo_channels_subset.h"

int main (int argc, char **argv) {
  /*
    Main program for generating the color combination files.
    Calls: alloc0 to allocate state,
           read_params to read input files,
	   build_coord_filenames to generate filenmes for its
      	     output files, the coordinate dominant combinations files,
              
	   alloc1 to allocate the subfields of state,
	   read_data to read the source, channels, and huebins input
	        files which should all be csv (comma-separated-values) files
	   prep to set up arrays and vectors that are fixed for 
	        the duration of execution.
           init_loop_state to sort and categorize the tsv values 
	        generated from the channels and source matrix.
	   gen_z_dom_combos to generate all possible combinations 
	        of z-dominant channels.
	   gen_y_dom_combos to generate all possible combinations 
	        of y_dominiant channels.
	   gen_x_dom_combos to generate all possible combinations 
	        of x-dominant channels.
  */
  struct state_struct *state;
  struct loop_state_struct *loop_state;
  double *curr_xyz;
  float  *tsv;
  char *parameter_file;      /* 128 */
  char *log_file;            /* 128 */
  char *output_file;         /* 128 */
  char *suffix_pos;          

  double test_rad;
  double u_target;
  double v_target;

  int num_channels;
  int nproc;

  int myid;
  int ierr;

  int success;
  int err_code;

  int fnl;
  int tsv_offset;

  int mask;
  int sbpi;

  FILE *lfp;
  FILE *ofp;

  lfp     = NULL;
  success = 1;
  nproc = 1;
  myid = 0;
  success = 1;
  /*
    Allocate a state structure.
  */
  success = alloc0(&state);
  if (success) {
    /*
      All procs need to open a log file.
    */
    log_file       = state->log_file;
    state->num_procs = nproc;
    state->my_id     = myid;
    parameter_file = state->parameter_file;
    /*
      Rank 0 reads in the parameters
    */
    if (argc < 2) {
      strcpy(parameter_file,"rthous_opt.in");
    } else {
      strcpy(parameter_file,argv[1]);
    }
    success = read_params(state);
  }
  if (success) {
#ifndef TM30
    /*
      If in select only mode (TM30 is not defined) then
      set select parameters that make sense.
    */
    state->select_only = 1;
#endif
    lfp = fopen(log_file,"w");
    state->lfp = lfp;
    success = alloc1(state);
  } else {
    fprintf(stderr,"gen_combos, read_params failed\n");
  }
  if (success) {
    num_channels = state->num_channels;
    /*
      May want to move the fill_isteps_per_channel call into init_loop_state
    */
    fill_isteps_per_channel(state);
    success = read_data(state);
    if (success == 0) {
      if (lfp) {
	fprintf(lfp,"gen_combos, Error reading data\n");
	fflush(lfp);
      }
    }
  } else {
    if (lfp) {
      fprintf(lfp,"gen_combos: alloc1 failed\n");
      fflush(lfp);
    }
  }
  /*
    Initalize polyshape info and 
    transpose the source and channels matrix as they are accessed
    column at a time.
    Form the sbar_t_channel and sbar_10_tchannel matrices,
    Form the mcat02 matrix and its lu factorization,
    Form the mhpe matrix,
    Set the scalar constants of state,
    Form the mu_wavelength, tplanck_num and c2_by_wavelength vectors,
    Form the table_planck matrix,
    Fill the chroma_shift_rr matrix
      
    NB. For rthous_opt we really only need to for sbar_t_channels matrix
    which is 3 by num_channels.
  */
  if (success) {
    build_coord_filenames(state);
    success = prep(state);
    /* check on isteps_per_channel vector, propogate last value
       for remaining channels if all values not set.
    */
    if (success == 0) {
      if (lfp) {
	fprintf(lfp,"gen_combos, error in prep\n");
	fflush(lfp);
      }
    }
  }
  if (success) {
    /* check for a subset, num_channels < full_num_channels
       If that is the case echo the subset channels to the subset_file
       (default is subset.csv)
    */
    if (state->num_channels < state->full_num_channels) {
      success = echo_channels_subset(state);
    }
  }      
  if (success) {
    /* 
       Set up the Nested loop over successive stages. 
       We loop in decreasing order.
    */
    if (lfp) {
      fprintf(lfp,"state->usage = %ld\n",state->usage);
    }
    num_channels = state->num_channels;
    u_target = state->u_target;
    v_target = state->v_target;
    test_rad = state->test_rad;
    success = init_loop_state(state,num_channels,u_target,v_target,
			      test_rad);
  }
  if (success) {
    /*
    success = gen_blue_combos(state);
    */
    success = bf_gen_z_dom_combos(state);
  }
  if (success) {
    /*
    success = gen_green_combos(state);
    */
    success = gen_y_dom_combos(state);
  }
  if (success) {
    /*
    success = gen_red_combos(state);
    */
    success = gen_x_dom_combos(state);
  }
  exit(1-success);
}
	  
      


