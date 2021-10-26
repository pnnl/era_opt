#ifndef _IVTBXR_H_
#define _IVTBXR_H_ 1
extern int ivtbxr(int n, int64_t *iv, int64_t *left_buffer, int64_t *right_buffer, 
		  int itag, struct comm_map_struct *reduction_map, FILE *lfp);
#endif
