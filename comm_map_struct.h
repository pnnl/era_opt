#ifndef _COMM_MAP_STRUCT_H_
#define _COMM_MAP_STRUCT_H_ 1
struct comm_map_struct {
  int parent;
  int left_child;
  int right_child;
  int my_rank;
  int nproc;
  int ierr;
  MPI_Request lc_event;
  MPI_Request rc_event;
  MPI_Request parent_event;
  MPI_Request quit_event;
}
;
#endif
