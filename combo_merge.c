#include "system_includes.h"
#include "combo_merge.h"
void  combo_merge(double *list1, 
		  double *list2, 
		  double *list3,
		  int l1,int l2, int nfields, int ifield) {
  /*
    Merge sorted combo records based on their
    ifield component  in ascending order.
    Called by: ts_sort
    Calls:     memcpy
  */
  int64_t move_size;
  int64_t imove_size;
  int64_t rec_size;
  int64_t l_8;
  double  *dp1;
  double  *dp2;
  double  *dp3;
  int j1;
  int j2;
  int j3;
  int n;
  l_8 = (int64_t)8;
  rec_size = nfields * l_8;
  n   = l1 + l2;
  j1 = 0;
  j2 = 0;
  j3 = 0;
  dp1 = list1;
  dp2 = list2;
  dp3 = list3;
  for (j3 = 0;j3 < n; j3++) {
    if (dp1[ifield] <= dp2[ifield]) {
      memcpy(dp3,dp1,rec_size);
      dp1 += nfields;       /* Caution address arithmetic */
      dp3 += nfields;       /* Caution address arithmetic */
      j1++;
      if (j1 == l1) {
	move_size  = (l2 - j2) * rec_size;
	if (move_size > 0) {
	  memcpy (dp3,dp2,move_size);
	}
	break;
      }
    } else {
      memcpy(dp3,dp2,rec_size);
      dp2 += nfields;      /* Caution address arithmetic */
      dp3 += nfields;      /* Caution address arithmetic */
      j2++;
      if (j2 == l2) {
	move_size = (l1-j1) * rec_size;
	if (move_size > 0) {
	  memcpy(dp3,dp1,move_size);
	}
	break;
      }
    }
  }
}
