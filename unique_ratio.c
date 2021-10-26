#include "system_includes.h"
#include "unique_ratio.h"
int unique_ratio(int num_channels, 
		 int max_index, 
		 int *istep_index,
		 int *unique_prime_factors, 
		 int *upf_index, 
		 int *mm_istep_index) {
  /*
    Determine if we have seen this ratio of steps before,
    if not set it to be its maximum possible multiple in the [1,max_index] range,
    and we return a value of 1,
    otherwise we return a value of 0.
    Called by: eprocess_combos
  */
  int i;
  int j;

  int min_step;
  int max_step;

  int istep;
  int ifact;

  int iqstep;
  int imult;

  int unique;
  int ipad;
  /*
    First we determine whether or not we have seen this ratio before.
    If the least common denominator (lcd) of the nonzero elements of
    istep_index is > 1 we have seen this combination before so
    return 0.
    Else return 1 and set mm_istep_index to be the maximum multiple
    of istep_index that fits in the power levels [1:max_index].
  */
  unique = 1;
  /*
    NB.  We should still improve on the efficiency here, by setting max_step and min_step
         to the first nonzero step value (instead of 1 and max_index), and then go back
         to only testing for the max if it wasn't the min.
  */
  min_step = max_index;
  max_step = 1;
  for (i=0;i<num_channels;i++) {
    istep = istep_index[i];
    if (istep > 0) {
      if (min_step > istep) {
	min_step = istep;
      } 
      if (max_step < istep) {
	max_step = istep;
      }
    }
  }
  if (min_step > 1) {
    /*
      This is what we term the lcd test (least common denominator).
      We need to test for divisibility by all the prime factors of 
      min_step.
    */
    for (j=upf_index[min_step];((j<upf_index[min_step+1]) && unique);j++) {
      unique = 0;
      ifact = unique_prime_factors[j];
      for (i=0;((i<num_channels) && (unique == 0));i++) {
	istep = istep_index[i];
	if (istep > min_step) {
	  iqstep = istep/ifact;
	  if ((iqstep*ifact) != istep) {
	    unique = 1;
	  }
	}
      }
    }
  }
  if (unique) {
    /* 
       Set mm_istep_index to be the maximum multple if istep_index in [1,max_index].
    */
    imult = max_index/max_step;
    if (imult == 1) {
      for (i=0;i<num_channels;i++) {
	mm_istep_index[i] = istep_index[i];
      }
    } else {
      for (i=0;i<num_channels;i++) {
	mm_istep_index[i] = istep_index[i] * imult;
      }
    }
  }
  return(unique);
}
