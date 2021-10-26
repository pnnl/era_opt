#include "system_includes.h"
#include "mpi.h"
#include "era_mpi_init.h"
int era_mpi_init(int *argc_p, char ***argv_p, int *myid_p, int *nproc_p) {
  /*
    Initialize parallel execution and set myid (rank) and nproc.
    Also set up MPI Error handler to return instead of abort
    so that relevant error messages may be printed.
    Called by: era_opt main program
    Calls:     MPI_Init, MPI_set_errhandler, MPI_Comm_rank,
               MPI_Comm_size, MPI_Abort, fprintf, fflush
  */
  int ierr;
  int success;
  int err_code;
  int myid;
  int nproc;
  int ipad;
  MPI_Comm comm;
  MPI_Comm ecomm;
  ierr = MPI_Init(argc_p,argv_p);
  comm = MPI_COMM_WORLD;
  /*
    Check on success of MPI_Init
  */
  success = 1;
  if (ierr != MPI_SUCCESS) {
    success = 0;
    fprintf(stderr,"rthous_opt: Error return from MPI_Init = %d\n",ierr);
    fflush(stderr);
    err_code   = 43;
    ierr = MPI_Abort(MPI_COMM_WORLD,err_code);
  }
  if (success) {
    /*
      Set mpi errors to return for testing and diagnostics instead of aborting:
      Always test return codes from MPI calls.
    */
    ierr = MPI_Comm_set_errhandler(comm,MPI_ERRORS_RETURN);
    if (ierr != MPI_SUCCESS) {
      success = 0;
      fprintf(stderr,"rthous_opt: Error return from "
	      "MPI_Comm_set_errhandler was %d\n",ierr);
      fflush(stderr);
      err_code   = 44;
      ierr = MPI_Abort(MPI_COMM_WORLD,err_code);
    }
  }
  /*
    Get local rank (myid) and number of processors (nproc)
  */
  if (success) {
    ierr = MPI_Comm_rank(comm,&myid);
    if (ierr != 0) {
      success = 0;
      fprintf(stderr,"rthous_opt: Error return from MPI_Comm_rank was %d\n",
	      ierr);
      fflush(stderr);
      err_code   = 45;
      ierr = MPI_Abort(MPI_COMM_WORLD,err_code);
    }
  }
  if (success) {
    ierr = MPI_Comm_size(comm,&nproc);
    if (ierr != 0) {
      success = 0;
      fprintf(stderr,"rthous_opt: Error return from MPI_Comm_size was %d\n",
	      ierr);
      fflush(stderr);
      err_code   = 46;
      ierr = MPI_Abort(MPI_COMM_WORLD,err_code);
    }
  }
  if (success) {
    *nproc_p = nproc;
    *myid_p  = myid;
  }
  return(success);
}
  
