#include "system_includes.h"
#include "comm_map_struct.h"
#include "init_ivtbar.h"
int init_ivtbar(int my_rank, int nproc, struct comm_map_struct **reduction_map_p, FILE *lfp) {
  /*
    Allocate a reduction_map.
    Set up a binary tree communication map in reduction_map.
    
    T = type   i = int, d = double, s = struct, f = file
    M = mode   0 = scalar, 1 = vector, 2 = array, * = pointer
    F = flow   i = input, o = output, b = input & output, w = scratch

    Variable   			Class  Description
               			TMF

    my_rank    			i0i    mpi rank of the calling process

    nproc                       i0i    total number of mpi ranks in the job.

    reduction_map_p             s*o    address of the allocated comm_map_struct.
                                       The parent, left_child, right_child, my_rank,
				       and nproc fields of the struct are set.

    lfp                         f*i    log file pointer for error messages.
				       
  */
  struct comm_map_struct *reduction_map;
  struct comm_map_struct reduction_map_i;
  int64_t one_l;
  int64_t ask_for;
  int success;
  int parent;
  int left_child;
  int right_child;
  success = 1;
  one_l = (int64_t)1;
  ask_for = (int64_t)sizeof(reduction_map_i);
  reduction_map = (struct comm_map_struct *)calloc(one_l,ask_for);
  if (reduction_map == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"init_ivtbar: Unable to allocate %ld bytes for reduction_map structure\n",
	      ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    *reduction_map_p = reduction_map;
    if (my_rank == 0) {
      parent = -1;
    } else {
      parent = (my_rank - 1) >> 1;
    }
    left_child = my_rank + my_rank + 1;
    right_child = left_child + 1;
    if (left_child >= nproc) {
      left_child = -1;
    }
    if (right_child >= nproc) {
      right_child = -1;
    }
    reduction_map->parent      = parent;
    reduction_map->left_child  = left_child;
    reduction_map->right_child = right_child;
    reduction_map->my_rank     = my_rank;
    reduction_map->nproc       = nproc;
  }
  return(success);
}
      
      
