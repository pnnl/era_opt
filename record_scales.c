#include "system_includes.h"
#include "pos_comma.h"
#include "record_scales.h"
void record_scales(char *value,int max_scales, int *number_scales_p,
		   double *scales, FILE *lfp) {
  /*
    Read a possibly commas separated list of values, might be just
    a single integer value.
    Called by: read_params
    Calls:     pos_comma,sscanf,fprintf,fflush
  */
  char *string;
  double scale;
  double prev_scale;

  int number_scales;
  int ns;

  int loc_comma;
  int ipad;

  string = value;
  number_scales = 0;
  prev_scale = 1.0;
  ns = sscanf(string,"%le",&scale);
  /*
  if (lfp) {
    fprintf(lfp,"record_scales: ns = %d\n",
	    ns);
    fflush(lfp);
  }
  */
  while((ns > 0) && (number_scales < max_scales)) {
    if ((scale >= 1.0) || (scale <= 0.0)) {
      if (lfp) {
	fprintf(lfp,"record_scales: Scales must be < 1.0 and > 0.0\n");
	fflush(lfp);
	number_scales = 0;
	ns = 0;
      }
    } else {
      if (scale >= prev_scale) {
	if (lfp) {
	  fprintf(lfp,"record_scales: Scales must strictly "
		  "decreasing in size\n");
	  fflush(lfp);
	  number_scales = 0;
	  ns = 0;
	}
      } else {
	scales[number_scales] = scale;
	number_scales += 1;
	ns = 0;
	loc_comma = pos_comma(string);
	string = string + loc_comma;
	if (loc_comma > 0) {
	  ns = sscanf(string,"%le",&scale);
	} else {
	  ns = 0;
	}
      }
    }
  }
  *number_scales_p  = number_scales;
}
