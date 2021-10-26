#include "system_includes.h"
#include "sort_combo_pairs_in_place.h"
#include "sort_combo_pairs.h"
#include "combo_merge.h"
#include "combo_sort.h"
void  combo_sort(int n, int nfields, int ifield, double *combos, 
		 double *scratch) {
  /*
    Sort a list of nx triples in xyz and their indices in imap
    by x if ifield = 0, by y if ifield = 1, and by z if ifield = 2.
    xyz_scratch and imap_scratch are scratch areas.
    The sorted list replaces the input.
    Called by: rthous_opt
    Calls:     sort_ts_pairs,sot_ts_pairs_in_place, ts_merge, memcpy
  */
  int64_t move_size;
  int64_t rec_size;
  int64_t l_8;
  double  *list1;
  double  *list2;
  double  *dtemp;

  int l1;
  int l2;

  int step;
  int j;

  int next_step;
  int ln;

  int k;
  int in_place_pairwise;

  int j1;
  int j2;

  int k1;
  int k2;
  
  list1  	    = combos;
  list2  	    = scratch;
  l_8    	    = (int64_t)8;
  rec_size          = nfields * l_8;
  in_place_pairwise = 0;
  if (n > 1) {
    for (step = 1; step < n; step += step) {
      in_place_pairwise = 1 - in_place_pairwise;
    }
    if (in_place_pairwise) {
      sort_combo_pairs_in_place(n,nfields,ifield,combos,scratch);
      list1 = combos;
      list2 = scratch;
    } else {
      sort_combo_pairs(n,nfields,ifield,combos,scratch);
      list1 = scratch;
      list2 = combos;
    }
  }
  for (step = 2; step < n; step += step) {
    next_step = step + step;
    for (j=0;j<(n-step);j+= next_step) {
      l1 = step;
      l2 = n - j - step;
      if (l2 > step) {
	l2 = step;
      }
      j1 = nfields*j;
      j2 = nfields*(j+step);
      combo_merge((double*)&list1[j1],
		  (double*)&list1[j2],
		  (double*)&list2[j1],
		  l1,l2,nfields,ifield);
    }
    ln = n & (next_step - 1);
    if (ln <= step) {
      move_size = ln * rec_size;
      k1 = n - ln;
      k2 = k1 * nfields;
      if (ln > 0) {
	memcpy((void*)&list2[k2],(void*)&list1[k2],move_size);
      }
    }
    dtemp  = list1;
    list1  = list2;
    list2  = dtemp;
  } /* end for step */
}
