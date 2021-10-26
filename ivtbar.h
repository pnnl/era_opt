#ifndef _IVTBAR_H_
#define _IVTBAR_H_ 1
extern int ivtbar(int n, 
		  int64_t *iv, 
		  int64_t *left_buffer, 
		  int64_t *right_buffer, 
		  int itag, 
		  struct comm_map_struct *reduction_map, 
		  FILE *lfp);
#endif
