#include "system_includes.h"
#include "ts_merge.h"
#include "sort_ts_pairs.h"
#include "sort_ts_pairs_in_place.h"
#include "ts_sort.h"
void  ts_sort(int n, int ifield, double *xyz, int *imap,
	      double *scratch,int *iscratch) {
  /*
    Sort a list of nx triples in xyz and their indices in imap
    by x if ifield = 0, by y if ifield = 1, and by z if ifield = 2.
    xyz_scratch and imap_scratch are scratch areas.
    The sorted list replaces the input.
    Called by: gen_combos
    Calls:     sort_ts_pairs,sot_ts_pairs_in_place, ts_merge, memcpy
  */
  int64_t move_size;
  int64_t imove_size;
  int64_t l_8;
  double *list1;
  double *list2;
  double *dtemp;
  int    *ilist1;
  int    *ilist2;
  int    *itemp;

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
  
  list1 = xyz;
  list2 = scratch;
  ilist1 = imap;
  ilist2 = iscratch;
  l_8    = (int64_t)8;
  
  in_place_pairwise = 0;
  if (n > 1) {
    for (step = 1; step < n; step += step) {
      in_place_pairwise = 1 - in_place_pairwise;
    }
    if (in_place_pairwise) {
      sort_ts_pairs_in_place(n,ifield,xyz,imap);
      list1 = xyz;
      ilist1 = imap;
      list2 = scratch;
      ilist2 = iscratch;
    } else {
      sort_ts_pairs(n,ifield,xyz,imap,scratch,iscratch);
      list1 = scratch;
      ilist1 = iscratch;
      list2 = xyz;
      ilist2 = imap;
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
      j1 = 3*j;
      j2 = 3*(j+step);
      ts_merge((double*)&list1[j1],(int*)&ilist1[j],
	       (double*)&list1[j2],(int*)&ilist1[j+step],
	       (double*)&list2[j1],(int*)&ilist2[j],
	       l1,l2,ifield);
    }
    ln = n & (next_step - 1);
    if (ln <= step) {
      move_size = 3*ln*l_8;
      imove_size = ln * 4;
      k1 = n - ln;
      k2 = k1 * 3;
      if (ln > 0) {
	memcpy((void*)&list2[k2],(void*)&list1[k2],move_size);
	memcpy((void*)&ilist2[k1],(void*)&ilist1[k1],imove_size);
      }
    }
    dtemp  = list1;
    list1  = list2;
    list2  = dtemp;
    itemp  = ilist1;
    ilist1 = ilist2;
    ilist2 = itemp;
  } /* end for step */
}
