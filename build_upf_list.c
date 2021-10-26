#include "system_includes.h"
#include "build_upf_list.h"
void build_upf_list(int max_index,
		    int mpf,
		    int *upf_index,
		    int *unique_prime_factors) {
  /*
    Build a list of the unique prime factors of the
    integer from 2 to max_index inclusive.

    upf_index must be of length is max_index + 2 elements.
    elements 0 and 1 are unused,
    for elements i in [2:max_index]
    upf_index[i] points to the start of integer i's unique prime factor
    list in the unique_prime_factors vector.
    upf_index[max_index+1] is one beyond max_index's last unique prime factor.
    So integer i's prime factors are in 
    unique_prime_factors[upf_index[i]:upf_index[i+1]-1] for i in [2:max_index]. 
    This array is used to count the number of unique factors when 
    building unique_prime_factors

    mpf is the maximum number of prime factors that could occur for any
    integer in [2:max_index] and unique_prime_factors must be 
    of length mpf * max_index.

    Called by init_loop_state.
  */
  int i;
  int ihalf_max_index;

  int inc_pos;
  int upf_pos;

  int j;
  int k;
  /*
    Use upf_index as a prime factor counter and also indicator of 
    primeness. Initialize it to all zero's to start with.
  */
  for (i=0;i<max_index+2;i++) {
    upf_index[i] = 0;
  } 
  /*
    Now we build a sieve of Erostosthenes to determine prime numbers and
    list the unique primes of numers in 2:max_index
    We only build up to max_index/2 as integers greater than max_index/2
    Will allready have had their primes listed if they are composit or
    have a count of 0 if they are prime at the end of this loop.
  */
  ihalf_max_index = max_index >> 1;
  for (i=2;i<=ihalf_max_index;i++) {
    if (upf_index[i] == 0) {
      /*
	This index has not yet been marked so i is prime.
	set its count to 1 and enter it into its own list.
      */
      inc_pos = i * mpf;
      upf_pos = inc_pos;
      unique_prime_factors[upf_pos] = i;
      upf_index[i] = 1;
      /*
	Now loop through and add this prime to all lists for all
	integer that are multiples of it, incrementing their prime
	factor count in the process.
      */
      for (j=i+i;j<=max_index;j+=i) {
	upf_pos += inc_pos;
	unique_prime_factors[upf_pos+upf_index[j]] = i;
	upf_index[j] += 1;
      }
    } /* end if i was prime upf_index[i] was 0. */
  }
  /* 
    Now go through and mark the primes that are greater than ihalf_max_index 
    Necessarily they won't have any multiples that are <= max_index.
  */
  upf_pos = ihalf_max_index * mpf;
  for (i=ihalf_max_index+1;i<=max_index;i++) {
    upf_pos += mpf;
    if (upf_index[i] == 0) {
      unique_prime_factors[upf_pos] = i;
      upf_index[i] = 1;
    }
  }
  /*
    Now upf_index[2:max_index] has the count of unique_primes for each
    integer in [2:max_index]. We partially sume the array giving pointers
    to where the unique prime factors for integer i are in the compressed
    array. We will comress the array after building this vector.
    j will be the partial sum.
  */
  j = 0;
  for (i=2;i<=max_index;i++) {
    k = upf_index[i];
    upf_index[i] = j;
    j = j + k;
  }
  upf_index[max_index+1] = j;
  /*
    Now compress the unique_prime_factors list in place.
    j will track the upf_index locations while upf_pos will
    track the ucompressed positions.
  */
  upf_pos = mpf;
  for (i=2;i<=max_index;i++) {
    upf_pos = upf_pos + mpf;
    k = 0;
    for (j=upf_index[i];j<upf_index[i+1];j++) {
      unique_prime_factors[j] = unique_prime_factors[upf_pos+k];
      k += 1;
    }
  }
}
