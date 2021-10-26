#include "system_includes.h"
#include "sort_ts_pairs.h"
void sort_ts_pairs(int n, int ifield, double *xyz, int *imap,
		   double *xyz2, int *imap2) {
  /*
    Sort pairs of tristimulus triples in descending order on ifield
    Called by: ts_sort
    Calls:     memcpy
  */
  double *dp1;
  double *dp2;
  double *dp3;
  int64_t l_24;
  int i;
  int nm1;
  int n_odd;
  int padi;
  l_24 = (int64_t)24;
  nm1 = n - 1;
  dp1 = xyz;
  dp2 = xyz + 3; /* Caution address arithmetic */
  dp3 = xyz2;
  for (i=0;i<nm1;i+=2) {
    if (dp1[ifield] < dp2[ifield]) {
      imap2[i] = imap[i+1];
      memcpy(dp3,dp2,l_24);
      dp3 += 3; /* Caution address arithmetic */
      imap2[i+1] = imap[i];
      memcpy(dp3,dp1,l_24);
      dp3 += 3; /* Caution address arithmetic */
    } else {
      imap2[i] = imap[i];
      memcpy(dp3,dp1,l_24);
      dp3 += 3; /* Caution address arithmetic */
      imap2[i+1] = imap[i+1];
      memcpy(dp3,dp2,l_24);
      dp3 += 3; /* Caution address arithmetic */
    }
    dp1 += 6;
    dp2 += 6;
  }
  n_odd = n & 1;
  if (n_odd) {
    imap2[n-1] = imap[n-1];
    memcpy(dp3,dp1,l_24);
  }
}
