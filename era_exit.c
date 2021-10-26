#include "system_includes.h"
#include "mpi.h"
#include "era_exit.h"
void era_exit(int nproc, int err_code) {
  int ierr;
  if (nproc > 1) {
    ierr = MPI_Abort(MPI_COMM_WORLD,err_code);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr,"MPI_Abort failed, err_code argument was %d\n",err_code);
      fflush(stderr);
    }
  } else {
    exit(err_code);
  }
}
