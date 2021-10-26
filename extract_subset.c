#include "system_includes.h"
#include "extract_subset.h"
void extract_subset(int num_channels_rows,int subset_size,int *subset_indices,
		    double *channels,double *channel_subset) {
  /*
    extract a subset of the columns of the channels matrix to use as 
    SPD for the thous_opt anlaysis - to allow mixing and matching parts
    of a large number of SPD's
    
    Called by: prep

    channels an array of SPD's stored in column major order.

    num_channels_rows is the length of a row in channels (and therefore also
    in channel_subset)
     
    subset_size is the number of columns to be extracted.
    subset_indices specify which columns (starting at 1) are to be extracted.
    
    channel_subset a channels_rows x subset_size matrix stored in column major order.

  */
  double *channels_row;
  double *subset_row;
  int k;
  int j;
  int i;
  subset_row = channel_subset;
  for (k=0;k<subset_size;k++) {
    j = (subset_indices[k] - 1) * num_channels_rows;
    channels_row = &channels[j];
    for (i=0;i<num_channels_rows;i++) {
      subset_row[i] = channels_row[i];
    }
    subset_row += num_channels_rows; /* Caution address arithmetic */
  }
}
