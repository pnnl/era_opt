#include "system_includes.h"
#include "state_struct.h"
#include "echo_channels_subset.h"
int echo_channels_subset(struct state_struct *state) {
  char *subset_file;
  double *channels;
  double *wavelength;
  double *channel_row;
  int num_channels;
  int channels_rows;
  int i;
  int j;
  int success;
  int ipad;
  FILE *lfp;
  FILE *sfp;
  success     	= 1;
  channels    	= state->channels;
  wavelength 	= state->wavelength;
  subset_file 	= state->subset_file;
  channels_rows = state->channels_rows;
  num_channels  = state->num_channels;
  lfp           = state->lfp;
  sfp = fopen(subset_file,"w+");
  if (sfp == NULL) {
    success =0;
    if (lfp) {
      fprintf(lfp,"echo_chanels_subset: Error, unable to open %s\n",
	      subset_file);
      fflush(lfp);
    }
  }
  if (success) {
    for(i=0;i<channels_rows;i++) {
      fprintf(sfp,"%lf",wavelength[i]);
      channel_row = &channels[i];
      for (j=0;j<num_channels;j++) {
	fprintf(sfp,",%lf",channel_row[0]);
	channel_row += channels_rows; /* Caution address arithmetic */
      }
      fprintf(sfp,"\n");
    }
    fclose(sfp);
  }
  return(success);
}
