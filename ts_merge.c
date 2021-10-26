#include "system_includes.h"
#include "ts_merge.h"
void  ts_merge(double *list1, int *ilist1,
	       double *list2, int *ilist2,
	       double *list3, int *ilist3,
	       int l1,int l2, int ifield) {
  /*
    Merge sorted tri-stimulus triples lists based on their
    ifield component (ifield = 0,1, or 2) in descending order.
    Called by: ts_sort
    Calls:     memcpy
  */
  int64_t move_size;
  int64_t imove_size;
  int64_t triple_size;
  int64_t l_4;
  int64_t l_8;
  double *dp1;
  double *dp2;
  double *dp3;
  int *ip1;
  int *ip2;
  int *ip3;
  int j1;
  int j2;
  int j3;
  int n;
  l_4 = (int64_t)4;
  l_8 = (int64_t)8;
  triple_size = 3 * l_8;
  n   = l1 + l2;
  j1 = 0;
  j2 = 0;
  j3 = 0;
  dp1 = list1;
  dp2 = list2;
  dp3 = list3;
  ip1 = ilist1;
  ip2 = ilist2;
  ip3 = ilist3;
  for (j3 = 0;j3 < n; j3++) {
    if (dp1[ifield] >= dp2[ifield]) {
      memcpy(dp3,dp1,triple_size);
      *ip3 = *ip1;
      dp1 += 3;       /* Caution address arithmetic */
      dp3 += 3;       /* Caution address arithmetic */
      ip1 += 1;       /* Caution address arithmetic */
      ip3 += 1;       /* Caution address arithmetic */
      j1++;
      if (j1 == l1) {
	imove_size = (l2 - j2) * l_4;
	move_size  = imove_size + imove_size;
	if (move_size > 0) {
	  memcpy (ip3,ip2,imove_size);
	  move_size = move_size * 3;
	  memcpy (dp3,dp2,move_size);
	}
	break;
      }
    } else {
      memcpy(dp3,dp2,triple_size);
      *ip3 = *ip2;
      dp2 += 3;      /* Caution address arithmetic */
      dp3 += 3;      /* Caution address arithmetic */
      ip2 += 1;      /* Caution address arithmetic */
      ip3 += 1;      /* Caution address arithmetic */
      j2++;
      if (j2 == l2) {
	imove_size = (l1-j1) * l_4;
	move_size  = imove_size + imove_size;
	if (move_size > 0) {
	  memcpy(ip3,ip1,imove_size);
	  move_size = move_size * 3;
	  memcpy(dp3,dp1,move_size);
	}
	break;
      }
    }
  }
}
