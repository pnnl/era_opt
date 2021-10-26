#include "system_includes.h"
#include "sort_ts_pairs_in_place.h"
void sort_ts_pairs_in_place(int n, int ifield,double *xyz, int* imap) {
  /*
    Called by: ts_sort
    Calls:     memcpy
  */
  double temp[3];
  double *dp1;
  double *dp2;
  int64_t l_24;
  int itemp;
  int i;
  int nm1;
  int padi;
  nm1 = n - 1;
  dp1 = xyz;
  dp2 = xyz + 3; /* Caution address arithmetic */
  l_24 = (int64_t)24;
  for (i=0;i<nm1;i+=2) {
    if (dp1[ifield] < dp2[ifield]) {
      itemp = imap[i];
      imap[i] = imap[i+1];
      imap[i+1] = itemp;
      memcpy(temp,dp1,l_24);
      memcpy(dp1,dp2,l_24);
      memcpy(dp2,temp,l_24);
    }
    dp1 += 6; /* Caution address arithmetic */
    dp2 += 6; /* Caution address arithmetic */
  }
}
	
