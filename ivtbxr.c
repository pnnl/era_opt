#include "system_includes.h"
#include "comm_map_struct.h"
#include "mpi.h"
#include "ivtbxr.h"
int ivtbxr(int n, int64_t *iv, int64_t *left_buffer, int64_t *right_buffer, 
	   int itag, struct comm_map_struct *reduction_map, FILE *lfp) {
  /*
    Integer Vector Tree-Based Maximum Reduction
    for an integer vector of length n, using provided buffer space
    and reduction communication map, we use a binary tree
    ipad is an alignment field. Parent = -1 means this node is the
    root of the communication tree. child = -1 is a nonexistant child.

    Called by: group_summary, lb_summary
    Calls:     MPI_Irecv, MPI_Wait, MPI_Send, fprintf, fflush

    T = type   i = int, d = double, s = struct, f = file
    M = mode   0 = scalar, 1 = vector, 2 = array, * = pointer
    F = flow   i = input, o = output, b = input & output, w = scratch

    Variable   			Class  Description
               			TMF

    n                           i0i    length in elements of the iv vector.

    iv                          i1b    local vector of vectors to be reduced.
                                       On rank 0 on ouput it has the final reduction.
				       On other ranks it has intermediate partial sums,
				       and should not be used.

    left_buffer                 i1w    buffer for receiving vector from left child.
                                       length is n.

    right_buffer                i1w    buffer for receiving vector from right child.

    itag                        i0i    message tag to be used in the point to point communications.
    
    reduction_map               s*i    communication map struct identifying parent, left-child and
                                       right-child fields.

    lfp                         f*i    log file pointer for error messages.				     
  */
  int left_child;
  int right_child;

  int parent;
  int i;

  int ierr;
  int success;

  int msg_type;
  int err_len;
  char *mpi_string;
  char err_buffer[512];
  
  MPI_Request left_event;
  MPI_Request right_event;
  MPI_Status  left_status;
  MPI_Status  right_status;
  MPI_Comm    comm;

  
  parent      = reduction_map->parent;
  left_child  = reduction_map->left_child;
  right_child = reduction_map->right_child;
  msg_type    = MPI_INT64_T;
  success     = 1;
  comm        = MPI_COMM_WORLD;
  mpi_string  = (char*)&err_buffer[0];
  /*
    Issue irecv to left child if it exists
  */
  if (left_child >= 0) {
    ierr = MPI_Irecv((void*)left_buffer, n, msg_type, left_child, itag, comm, &left_event);
    if (ierr != MPI_SUCCESS) {
      success = 0;
      if (lfp) {
	MPI_Error_string(ierr,mpi_string,&err_len);
	fprintf(lfp,"ivtbxr: Error return from MPI_Irecv from %d was %d, %s\n",left_child,
		ierr, mpi_string);
	fflush(lfp);
      }
    }
  }
  /*
    Issue irecv to right child if it exists
  */
  if (success) {
    if (right_child >= 0) {
      ierr = MPI_Irecv((void*)right_buffer, n, msg_type, right_child, itag, comm, &right_event);
      if (ierr != MPI_SUCCESS) {
	success = 0;
	if (lfp) {
	  MPI_Error_string(ierr,mpi_string,&err_len);
	  fprintf(lfp,"ivtbxr: Error return from MPI_Irecv from %d was %d, %s\n",right_child,
		  ierr, mpi_string);
	  fflush(lfp);
	}
      }
    }
  }
  if (success) {
    if (left_child >= 0) {
      ierr = MPI_Wait(&left_event,&left_status);
      if (ierr != MPI_SUCCESS) {
	success = 0;
	if (lfp) {
	  MPI_Error_string(ierr,mpi_string,&err_len);
	  fprintf(lfp,"ivtbxr: Error return from MPI_Wait from %d was %d, %s\n",left_child,ierr,
		  mpi_string);
	  fflush(lfp);
	}
      } else {
	/*
	  Combine the left child results
	*/
	for (i=0;i<n;i++) {
	  if (left_buffer[i] > iv[i]) {
	    iv[i] = left_buffer[i];
	  }
	}
      }
    } /* end if (left_child > 0) */
  } /* end if success */
  if (success) {
    if (right_child >= 0) {
      ierr = MPI_Wait(&right_event,&right_status);
      if (ierr != MPI_SUCCESS) {
	success = 0;
	if (lfp) {
	  MPI_Error_string(ierr,mpi_string,&err_len);
	  fprintf(lfp,"ivtbxr: Error return from MPI_Wait from %d was %d, %s\n",right_child,ierr,
		  mpi_string);
	  fflush(lfp);
	}
      } else {
	/*
	  Combine the right child results
	*/
	for (i=0;i<n;i++) {
	  if (right_buffer[i] > iv[i]) {
	    iv[i] = right_buffer[i];
	  }
	}
      }
    } /* end if (right_child > 0) */
  } /* end if success */
  if (success) {
    if (parent >= 0) {
      ierr = MPI_Send(iv,n,msg_type,parent,itag,comm);
      if (ierr != MPI_SUCCESS) {
	success = 0;
	if (lfp) {
	  MPI_Error_string(ierr,mpi_string,&err_len);
	  fprintf(lfp,"ivtbxr: Error return from MPI_Send to %d was %d,%s\n",parent,ierr,
		  mpi_string);
	  fflush(lfp);
	}
      }
    }
  }
  return(success);
}
