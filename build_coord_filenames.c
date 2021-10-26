#include "system_includes.h"
#include "state_struct.h"
#include "extract_root.h"
#include "add_suffix.h"
#include "build_coord_filenames.h"
void build_coord_filenames(struct state_struct *state) { 
  /*
    Build the z_dom_file, y_com_file, and x_dom_file file names from
    the root of the channels_file file name, if they are empty.
    Called by: gen_combos, era_opt
    Calls:     extract_root,strcpy,add_suffix,fprintf, fflush
  */
  char *z_dom_file;
  char *y_dom_file;
  char *x_dom_file;
  char *summary_file;
  char *channels_file;
  int  extra_log_info;
  int  ipad;
  FILE *lfp;
  FILE *efp;
  z_dom_file    = state->z_dom_file;
  y_dom_file    = state->y_dom_file;
  x_dom_file    = state->x_dom_file;
  channels_file = state->channels_file;
  /*
  summary_file  = state->summary_file;
  */
  extra_log_info = state->extra_log_info;
  if (strlen(z_dom_file) == 0) {
    extract_root(channels_file,z_dom_file);
    add_suffix(z_dom_file,"z.dat");
  }
  if (strlen(y_dom_file) == 0) {
    extract_root(channels_file,y_dom_file);
    add_suffix(y_dom_file,"y.dat");
  }
  if (strlen(x_dom_file) == 0) {
    extract_root(channels_file,x_dom_file);
    add_suffix(x_dom_file,"x.dat");
  }
  /*
  extract_root(channels_file,summary_file);
  add_suffix(summary_file,"sum");
  */
  if (extra_log_info) {
    lfp = state->lfp;
    if (lfp) {
      fprintf(lfp,"z_dom_file   = %s\n",z_dom_file);
      fprintf(lfp,"y_dom_file   = %s\n",y_dom_file);
      fprintf(lfp,"x_dom_file   = %s\n",x_dom_file);
      /*
      fprintf(lfp,"summary_file = %s\n",summary_file);
      */
      fflush(lfp);
    }
  }
}
