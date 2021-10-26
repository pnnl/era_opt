#include "system_includes.h"
#include "state_struct.h"
#include "groups_struct.h"
#include "comm_map_struct.h"

#include "era_mpi_init.h"
#include "era_exit.h"
#include "alloc0.h"
#include "read_params.h"
#include "load_params.h"
#include "dist_params.h"
#include "alloc1.h"
#include "init_ivtbar.h"
#include "fill_isteps_per_channel.h"
#include "read_data.h"
#include "dist_data.h"
#include "build_coord_filenames.h"
#include "prep.h"
#include "eprocess_combos.h"
#include "bf_process_combos.h"
#ifdef TM30
#include "alloc3.h"
#include "init_label_maps.h"
#include "rprocess_combos.h"
#include "group_summary.h"
#include "print_summary.h"
#include "select_summary.h"
#include "cct_rf_summary.h"
#endif
#include "lb_summary.h"
/*
  Add a -DDBG="1" to compile flags instead of uncommenting this.
#define DBG 1
*/
#ifndef DBG
#include "mpi.h"
#endif
#include  "init_loop_state.h"
#include  "mmap_coord_files.h"
int main (int argc, char **argv) {
  /*
    Main program for the era_opt code that
    originated from Thous_opt_4.m matlab code from Alp Durmus.
    Calls: alloc0 to allocate state,
           read_params to read input files,
	   alloc1 to allocate the subfields of state,
	   read_data to read the source, channels, and huebins input
	        files which should all be csv (comma-separated-values) files
	   prep to set up arrays and vectors that are fixed for 
	        the duration of execution.
	   
	   ies_tm_30_18.h for combos that pass the 3500K Planckian test
	   
  */
  struct state_struct *state;
  struct groups_struct *groups_state;
  struct comm_map_struct *reduction_map;
  struct loop_state_struct *loop_state;
  char *parameter_file;      /* 128 */
  char *log_file;            /* 128 */
  char *output_file;         /* 128 */
  char *combos_out_file;     /* 128 */
  char *summary_file;        /* 128 */
  char *output_dir;          /* 128 */
  char *dir_plus_fn;         /* 256 */
  char *suffix_pos;          
  char *file_pos;
  int  *isteps_per_channel;

  double test_rad;
  double u_target;
  double v_target;
  double rparams[8];
  
  int iparams[64];

  int num_channels;
  int nproc;

  int myid;
  int ierr;

  int success;
  int err_code;

  int fnl;
  int lfn;

  int i;
  int num_huebins;

  int full_output;
  int select_only;

  int iout_fn_len;
  int iout_dir_len;

  int nf;
  int read_combos;

  int dump_combos;
  int ipad;

  FILE *lfp;
  FILE *ofp;
  FILE *sfp;
  FILE *dfp;
  lfp     = NULL;
  success = 1;
#ifdef DBG
  nproc = 1;
  myid = 0;
  success = 1;
#else 
  success = era_mpi_init(&argc,&argv,&myid,&nproc);
#endif
  /*
    All processors allocate a state structure, and
    space for the filenames.
  */
  if (success) {
    success = alloc0(&state);
  } else {
    fprintf(stderr,"era_mpi_init failed\n");
    err_code = 47;
    era_exit(nproc,err_code);
  }
  if (success) {
    /*
      All procs need to open a log file.
    */
    log_file       = state->log_file;
    state->num_procs = nproc;
    state->my_id     = myid;
    parameter_file = state->parameter_file;
    if (myid == 0) {
      /*
	Rank 0 reads in the parameters
      */
      if (argc < 2) {
	strcpy(parameter_file,"era_opt.in");
      } else {
	strcpy(parameter_file,argv[1]);
      }
      success = read_params(state);
#ifndef TM30
      /*
	If in select only mode (TM30 is not defined) then
	set select parameters that make sense.
      */
      state->select_only = 1;
      state->dump_combos = 1;
      state->full_output = 0;
#endif
      if (success) {
	build_coord_filenames(state);
      }
      /*
	load the scalar int and double parameters into the iparams
	and rparams vectors, include success as the first integer parameter.
      */
      load_params(success,state,iparams,rparams);
    }
  } else {
    iparams[0] = 0;
    success    = 0;
    err_code   = 48;
    fprintf(stderr,"era_opt, alloc0 failed\n");
    fflush(stderr);
    era_exit(nproc,err_code);
  }
  /*
    Now we need to broadcast out all the filenames, 
    iparams and rparams vectors,
    This happens in dist_params.
  */
  fnl       = 128; /* Filename lengths */
  nf        = 32;  /* Number of filenames/directories allocated in alloc0 */
  if (nproc > 1) {
    success = dist_params(state,fnl,nf,nproc,iparams,rparams);
  }
  full_output = state->full_output;
  select_only = state->select_only;
  read_combos = state->read_combos;
  dump_combos = state->dump_combos;
  /*
    Build and open a unique log file for each proc.
  */
  if (success) {
    lfn = strlen(log_file);
    if (lfn + 8 > fnl-1) {
      /* Leave at least 9 characters for suffix, 
	 Truncate if too long. */
      lfn = fnl-9;
      log_file[lfn] = '\0';
    }
    suffix_pos =  log_file + lfn; /* Caution addres arith.*/
    sprintf(suffix_pos,".%d",myid);
    
    lfp = fopen(log_file,"w");
    state->lfp = lfp;
    /*
      At this point the broadcast of the parameters has succeeded.
      Allocate space for the input data and setup data derivations
    */
    success = alloc1(state);
  } else {
    err_code = 52;
    fprintf(stderr,"era_opt, dist_params failed\n");
    fflush(stderr);
    era_exit(nproc,err_code);
  }
  if (success) {
    if (select_only == 0) {
      /*
	Allocate space for and initialize the groups_state structure.
      */
#ifdef TM30
      success = alloc3(state,&groups_state);
      state->usage += groups_state->usage;
      if (success) {
	state->groups_state = groups_state;
	init_label_maps(groups_state);
      }
#endif
    }
  }
  if (success) {
    /*
      Allocate space for and initialize the the reduction communication map structure.
    */
    success = init_ivtbar(myid,nproc,(struct comm_map_struct **)&reduction_map,lfp);
    if(success) {
      state->reduction_map = reduction_map;
    }
  }
  if (success) {
    num_channels = state->num_channels;
    num_huebins  = state->num_huebins;
    /*
      if (lfp) {
      fprintf(lfp,"nc_steps_set = %d\n",nc_steps_set);
      fprintf(lfp,"num_channels = %d\n",num_channels);
      fflush(lfp);
      }
    */
    if (myid == 0) {
      /*
	May want to move the fill_isteps_per_channel call into init_loop_state
      */
      fill_isteps_per_channel(state);
      success = read_data(state);
      if (success == 0) {
	if (lfp) {
	  fprintf(lfp,"era_opt, Error reading data on rank 0\n");
	  fflush(lfp);
	}
	err_code = 61;
	era_exit(nproc,err_code);
      }
    }
  } else {
    success = 0;
    if (lfp) {
      fprintf(lfp,"era_opt: alloc1 failed\n");
      fflush(lfp);
    }
    err_code = 61;
    era_exit(nproc,err_code);
  }
  /* 
    Now we need to broadcast the channel_matrix, isteps_per_channel,
    subset_indices, test_radii_sq, the source_matrix, 
    the huebins, and the cct_matrix, to all ranks.
    The latter 3 only if select_only == 0.
  */
  if (success) {
    if (nproc > 1) {
      success = dist_data(state);
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
      
    NB. For select_only we really only need to form sbar_t_channels matrix
    which is 3 by num_channels.
  */
  if (success) {
    success = prep(state);
    /* check on isteps_per_channel vector, propogate last value
       for remaining channels if all values not set.
    */
    if (success == 0) {
      if (lfp) {
	fprintf(lfp,"era_opt, error on rank 0 in prep\n");
	fflush(lfp);
      }
      err_code = 81;
      era_exit(nproc,err_code);
    }
  }
  if (success) {
    output_file  = state->output_file;
    iout_fn_len = strlen(output_file);
    suffix_pos =  output_file + iout_fn_len; /* Caution address arithmetic */
    sprintf(suffix_pos,".%d\0",myid);
    if (full_output) {
      /*
	Prepend output_dir to output_file 
      */
      output_dir  = state->output_dir;
      dir_plus_fn = state->dir_plus_fn;
      file_pos = dir_plus_fn;
      iout_dir_len = strlen(output_dir);
      if (iout_dir_len > 0) {
	strcpy(dir_plus_fn,output_dir);
	file_pos = dir_plus_fn + iout_dir_len; /* Caution address arithmetic */
      }
      strcpy(file_pos,output_file);
      ofp = fopen(dir_plus_fn,"w");
      if (ofp == NULL) {
	success = 0;
	fprintf(stderr,"era_opt Error: could not open output file %s\n",
		dir_plus_fn);
	fflush(stderr);
	err_code = 82;
	era_exit(nproc,err_code);
      } else {
	state->ofp = ofp;
      }
    }
  }
  if (success) {
    if (select_only && dump_combos) {
      combos_out_file = state->combos_out_file;
      iout_fn_len = strlen(combos_out_file);
      suffix_pos =  combos_out_file + iout_fn_len;
      sprintf(suffix_pos,".%d\0",myid);

      output_dir  = state->output_dir;
      dir_plus_fn = state->dir_plus_fn;
      file_pos = dir_plus_fn;
      iout_dir_len = strlen(output_dir);
      if (iout_dir_len > 0) {
	strcpy(dir_plus_fn,output_dir);
	file_pos = dir_plus_fn + iout_dir_len; /* Caution address arithmetic */
      }
      strcpy(file_pos,combos_out_file);
      dfp = fopen(dir_plus_fn,"w");
      if (dfp == NULL) {
	success = 0;
	fprintf(stderr,"era_opt Error: could not open combos_out_file %s\n",
		combos_out_file);
	fflush(stderr);
	err_code = 82;
	era_exit(nproc,err_code);
      } else {
	state->dfp = dfp;
      }
    }
  }
  if (success) {
    /*
      Only proc 0 opens a summary file.
    */
    if (myid == 0) {
      summary_file = state->summary_file;
      /*
	Prepend output_dir to summary file name .
      */
      output_dir  = state->output_dir;
      iout_dir_len = strlen(output_dir);
      dir_plus_fn = state->dir_plus_fn;
      file_pos = dir_plus_fn;
      if (iout_dir_len > 0) {
	strcpy(dir_plus_fn,output_dir);
	file_pos = dir_plus_fn + iout_dir_len; /* Caution address arithmetic */
      }
      strcpy(file_pos,summary_file);
      /*
      sfp = fopen(summary_file,"w");
      */
      sfp = fopen(dir_plus_fn,"w");
      if (sfp == NULL) {
	success = 0;
	fprintf(stderr,"era_opt Error: could not open summary file %s\n",
		summary_file);
	fflush(stderr);
	err_code = 921;
	era_exit(nproc,err_code);
      } else {
	state->sfp = sfp;
      }
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
  /*
  ierr=MPI_Barrier(MPI_COMM_WORLD);
  */
  if (success) {
    success = mmap_coord_files(state);
  }
  if (success) {
    loop_state = state->loop_state;
    loop_state->ofp = ofp;
    loop_state->lfp = lfp;
    if (full_output) {
      /*
	Print a label row in the log file.
	We no longer print all the output to the log file no need
	for a header.
      */
      /*
      for (i=1;i<=num_channels;i++) {
	fprintf(lfp,"S%2d,",i);
      }
      fprintf(lfp,"r_f,r_g,LER,M/P,SUM(P),RF_H1,PHI,SVD_PHI,ELL_STAT,");
      fprintf(lfp,"NINT_RF,NINT_RG,LGUD,L9UC,L9UD,L5UC,");
      fprintf(lfp,"L5UO,L9BD,L5BC,L5BD,E9BD,E5BC,E5BD,");
      fprintf(lfp,"P,V,F,S9BD,S5BC,S5BD,DIST,NZC\n");
      */
      /*
	Print a label row in the output file.
      */
      for (i=1;i<=num_channels;i++) {
	fprintf(ofp,"P%2d,",i);
      }
      fprintf(ofp,"r_f,r_g,");
      for (i=1;i<=num_huebins;i++) {
	fprintf(ofp,"RCS_H%2d,",i);
      }
      /*
      fprintf(ofp,"LER,M/P,SUM(P),RF_H1,PHI,SVD-PHI,ST,");
      */
      fprintf(ofp,"LER,M/P,SUM(P),RF_H1,PHI,ST,");
      fprintf(ofp,"NINT_RF,NINT_RG,LGUD,L9UC,L9UD,L5UC,");
      fprintf(ofp,"L5UO,L9BD,L5BC,L5BD,E9BD,E5BC,E5BD,");
      /*
      fprintf(ofp,"P,V,F,S9BD,S5BC,S5BD,DIST\n");
      */
      fprintf(ofp,"P,V,F,DIST,CCT_X,CCT_Y,CCT_Z,T_T\n");
    }
    if (select_only && dump_combos) {
      isteps_per_channel = state->isteps_per_channel;
      /*
	Here we are assuming all channels have the same number
	of power steps.
      */
      fprintf(dfp,"%d,%d\n",num_channels,isteps_per_channel[0]);
      fflush(dfp);
    }
    if (read_combos != 0) {
#ifdef TM30
#ifdef DBG
      rprocess_combos(state,loop_state,myid,nproc);
#endif
#endif
    } else {
      if (select_only == 0) {
	eprocess_combos(state,loop_state,myid,nproc);
      } else {
	if (state->use_bruteforce == 0) {
	  eprocess_combos(state,loop_state,myid,nproc);
	} else {
	  bf_process_combos(state,loop_state,myid,nproc);
	}
      }
    }
  }
  if (success) {
    if (select_only == 0) {
#ifdef TM30
      group_summary(nproc,state);
      cct_rf_summary(nproc,state);
      if (myid == 0) {
	print_summary(state,sfp);
      }
#endif
    }
    lb_summary(nproc,state);
    if (myid == 0) {
      fprintf(sfp," Load balance information:\n");
      fprintf(sfp," number of tasks, izc                       = %ld\n",state->izc);
      fprintf(sfp," number of subtasks iy_delta_sum            = %ld\n",state->iy_delta_sum);
      fprintf(sfp," number of combinations examined            = %ld\n",state->ix_delta_sum);
      fprintf(sfp," total combinations in TCA                  = %ld\n",state->num_in_tca);
      fprintf(sfp," max subtasks/node(max_iy_delta_sum         = %ld\n",state->max_iy_delta_sum);
      fprintf(sfp," max combos examined/node(max_ix_delta_sum) = %ld\n",state->max_ix_delta_sum);
      fprintf(sfp," max combos in TCA/node (max_in_tca)        = %ld\n",state->max_in_tca);
      if (select_only == 0) {
         fprintf(sfp," Total number of results                    = %ld\n",state->total_number_results);
         fprintf(sfp," max number results/node                    = %ld\n",state->max_number_results);
      }
      fflush(sfp);
      fclose(sfp);
    }
  } 
  if (lfp) {
    fclose(lfp);
  }
  if (full_output) {
    if (ofp) {
      fclose(ofp);
    }
  }
  if (select_only & dump_combos) {
    if (dfp) {
      fclose(dfp);
    }
  }
#ifndef DBG  
  if (nproc > 1) {
    ierr=MPI_Finalize();
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr,"Error return from MPI_Finalize was %d\n",ierr);
      fflush(stderr);
    }
  }
#endif
  return(0);
}
