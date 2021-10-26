#include "system_includes.h"
#include "pos_comma.h"
#include "record_subset.h"
void record_subset(char *value, int max_subset_size, int num_channels, 
		   int *subset_size_p, int *subset_indices, FILE *lfp) {
  /*
    Read a possibly commas separated list of values, might be just
    a single integer value.
    Called by: read_params
    Calls:     pos_comma
  */
  char *string;

  int subset_size;
  int nv;

  int ns;
  int loc_comma;

  string = value;
  subset_size = 0;
  ns = sscanf(string,"%d",&nv);
  /*
  if (lfp) {
    fprintf(lfp,"record_subset: ns = %d, nv = %d\n",
	    ns,nv);
    fflush(lfp);
  }
  */
  while((ns > 0) && (subset_size < max_subset_size)) {
    if (nv > 0) {
      if (nv <= num_channels) {
	subset_indices[subset_size] = nv;
	subset_size++;
      } else {
	fprintf(lfp,"record_subset: Error subset indices must be "
		"in [1,num_channels], subset index of %d found\n",nv);
	fflush(lfp);
      }
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
  *subset_size_p = subset_size;
}
