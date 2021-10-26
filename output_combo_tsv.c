#include "system_includes.h"
#include "output_combo_tsv.h"
void output_combo_tsv(FILE *dfp, int num_channels, 
		      int *ipstep_index, 
		      double x, double y, double z) {
  /*
    Print the integer power steps (unperumted) in ipstep_index
    and the tristiumulus (x,y,z) values for a combination to the
    combos_out_file.
  */
  int i;
  if (dfp) {
    for(i=0;i<num_channels;i++) {
      fprintf(dfp,"%3d,",ipstep_index[i]);
    }
    fprintf(dfp,"%le,%le,%le\n",x,y,z);
  }
  /* we do not flush dfp, to avoid many small writes */
}
