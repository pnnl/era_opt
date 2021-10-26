#include "system_includes.h"
#include "count_nzc.h"
int count_nzc(int num_channels, int *istep_index) {
  /*
    Count the number of nonzer power levels in istep_index vector
    of length num_channels.

    Called by: process_combos, ms_process_combos;
  */
  int nzc;
  int i;
  nzc = 0;
  for (i=0;i<num_channels;i++) {
    if (istep_index[i] > 0) nzc += 1;
  }
  return(nzc);
}
