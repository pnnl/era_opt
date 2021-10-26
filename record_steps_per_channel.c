#include "system_includes.h"
#include "pos_comma.h"
#include "record_steps_per_channel.h"
void record_steps_per_channel(char *value,int *nc_steps_set_p,
			      int *isteps_per_channel, FILE *lfp) {
  /*
    Read a possibly commas separated list of values, might be just
    a single integer value.
    Called by: read_params
    Calls:     pos_comma
  */
  char *string;

  int nc_steps_set;
  int nv;

  int ns;
  int loc_comma;

  string = value;
  nc_steps_set = 0;
  ns = sscanf(string,"%d",&nv);
  /*
  if (lfp) {
    fprintf(lfp,"record_steps_per_channel: ns = %d, nv = %d\n",
	    ns,nv);
    fflush(lfp);
  }
  */
  while(ns > 0) {
    if (nv > 0) {
      isteps_per_channel[nc_steps_set] = nv;
      nc_steps_set++;
      ns = 0;
      loc_comma = pos_comma(string);
      string = string + loc_comma;
      if (loc_comma > 0) {
	ns = sscanf(string,"%d",&nv);
      } else {
	ns = 0;
      }
    } else {
      ns = 0;
    }
  }
  *nc_steps_set_p = nc_steps_set;
}
