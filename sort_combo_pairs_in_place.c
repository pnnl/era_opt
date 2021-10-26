#include "system_includes.h"
#include "sort_combo_pairs_in_place.h"
void sort_combo_pairs_in_place(int n, int nfields, 
			       int ifield, double *combos, double *scratch) {
  /*
    Sort combination pairs in increasing order in ifield
    Called by: combo_sort
    Calls:     memcpy
  */
  double *dp1;
  double *dp2;
  int64_t rec_len;
  int itemp;
  int i;
  int nm1;
  int two_nfields;
  nm1 = n - 1;
  dp1 = combos;
  dp2 = combos + nfields; /* Caution address arithmetic */
  two_nfields = nfields + nfields;
  rec_len = (int64_t)(nfields * 8);
  for (i=0;i<nm1;i+=2) {
    if (dp1[ifield] > dp2[ifield]) {
      memcpy(scratch,dp1,rec_len);
      memcpy(dp1,dp2,rec_len);
      memcpy(dp2,scratch,rec_len);
    }
    dp1 += two_nfields; /* Caution address arithmetic */
    dp2 += two_nfields; /* Caution address arithmetic */
  }
}
	
