#include "system_includes.h"
#include "pos_comma.h"
#include "record_dbg_steps.h"
void record_dbg_steps(char *value,int num_channels, int *num_dbg_steps,
		      double *dbg_steps, FILE *lfp) {
  /*
    Read a comma separated list of values
    Called by: read_params
    Calls:     pos_comma,sscanf,fprintf,fflush
  */
  char *string;
  double step;
  double prev_scale;

  int ndbgs;
  int ns;

  int loc_comma;
  int i;
  

  ndbgs = 0;
  string = value;
  for (i=0;i<num_channels;i++) {
    dbg_steps[i] = 0.0;
  }
  for (i=0;i<num_channels;i++) {
    ns = sscanf(string,"%le",&step);
    if (ns == 1) {
      ndbgs += 1;
      dbg_steps[i] = step;
    } else {
      /*
	do we want to print a warning to lfp here?
      */
      dbg_steps[i] = 0.0;
    }
    loc_comma = pos_comma(string);
    string = string + loc_comma;
    if (loc_comma == 0) break;
  }
  *num_dbg_steps = ndbgs;
}
