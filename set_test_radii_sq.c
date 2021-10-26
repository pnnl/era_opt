#include "system_includes.h"
#include "set_test_radii_sq.h"
int set_test_radii_sq(int num_channels, int use_vrb, int num_vrb,
		      int *vrb_first_index, int *vrb_last_index,
		      double *vrb_radii, double test_rad,
		      double *test_radii_sq) {
  /*
    Determine the square of the test_radii for each of the different
    number of nonzero channels from the vrb_first_index,
    vrb_last_index, and vrb_radii vectors. But if use_vrb is 0,
    then set all values [1:num_channels] of test_radii_sq to be
    test_rad * test_rad;

    Arguments         TMF     Description
    num_channels      I0I     number of channels
    
    use_vrb           I0I     0 for fixed test radius, test_rad
                              1 for variable test radius based
			      on number of nonzero channels.

			      
    num_vrb           I0I     number of variable radius bins
                              vrb_first_index, vrb_last_index
			      and vrb_radii are expected to
			      have elements [1:num_channels] set
			      if use_vrb is 1.

    vrb_first_index   I1I     first index of variable bins (inclusive)
                              Length is num_channels + 1 ( as we
			      use elements [1:num_vrb]

    vrb_last_index    I1I     last index of variable bins (inclusive)
                              Length is num_channels + 1 ( as we
			      use elements [1:num_vrb]

    vrb_radii         D1I     test radii to use for numbers of nonaero
                              channels delimiet by the vrb_first_index
			      and vrb_last_index arrays.

    test_rad          D0I     Test radius to be used for all numbers
                              of channels if use_vrb is 0.

    test_radii_sq     D1O     vector of length 1 + num_channels whose
                              elements [1:num_channels] is set to
			      the square of the test_radius for 
			      combinations with the index number of 
			      nonzero power levels for channels.

    Called by: read_params.			      

  */
  double test_rad_sq;
  int i;
  int j;
  int success;
  int ipad;
  success = 1;
  if (use_vrb == 0) {
    test_rad_sq = test_rad * test_rad;
    for (i=1;i<=num_channels;i++) {
      test_radii_sq[i] = test_rad_sq;
    }
  } else {
    /*
      Verify that (vrb_first_index[i],vrb_last_index[i]) for i=1,
      num_vrb 
    */
    if (vrb_first_index[0] != 1) {
      success = 0;
    } else {
      for (i=0;((i<num_vrb-1) && success);i++) {
	if ((vrb_last_index[i]+1) != vrb_first_index[i+1]) {
	  success = 0;
	}
      }
      if (vrb_last_index[num_vrb-1] < num_channels) {
	success = 0;
      }
    }
    if (success) {
      for (i=0;i<num_vrb;i++) {
	test_rad_sq = vrb_radii[i];
	test_rad_sq = test_rad_sq * test_rad_sq;
	for (j=vrb_first_index[i];j<=vrb_last_index[i];j++) {
	  test_radii_sq[j] = test_rad_sq;
	}
      }
    }
  }
  return(success);
}
