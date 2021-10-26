#include "system_includes.h"
#include "sort_combo_pairs.h"
void sort_combo_pairs(int n, int nfields, int ifield, double *combos, 
		      double *combos2) {
  /*
    Sort combinations in increasing order on ifield
    Called by: combo_sort
    Calls:     memcpy
  */
  double  *dp1;
  double  *dp2;
  double  *dp3;
  int64_t rec_len;
  int i;
  int nm1;
  int n_odd;
  int two_nfields;
  rec_len = (int64_t)(nfields*8);
  nm1 = n - 1;
  two_nfields = nfields+nfields;
  dp1 = combos;
  dp2 = combos + nfields; /* Caution address arithmetic */
  dp3 = combos2;
  for (i=0;i<nm1;i+=2) {
    if (dp1[ifield] > dp2[ifield]) {
      memcpy(dp3,dp2,rec_len);
      dp3 += nfields; /* Caution address arithmetic */
      memcpy(dp3,dp1,rec_len);
      dp3 += nfields; /* Caution address arithmetic */
    } else {
      memcpy(dp3,dp1,rec_len);
      dp3 += nfields; /* Caution address arithmetic */
      memcpy(dp3,dp2,rec_len);
      dp3 += nfields; /* Caution address arithmetic */
    }
    dp1 += two_nfields;
    dp2 += two_nfields;
  }
  n_odd = n & 1;
  if (n_odd) {
    memcpy(dp3,dp1,rec_len);
  }
}
