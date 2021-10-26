#ifndef _DVTBNR_H_
#define _DVTBNR_H_ 1
extern int dvtbnr(int n, double *dv, double *left_buffer, double *right_buffer, 
		  int itag, struct comm_map_struct *reduction_map, FILE *lfp);
#endif
