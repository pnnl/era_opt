#ifndef _RECORD_SUBSET_H_
#define _RECORD_SUBSET_H_ 1
extern void record_subset(char *value, int max_subset_size, int num_channels, 
			  int *subset_size_p, int *subset_indices, FILE *lfp);
#endif
