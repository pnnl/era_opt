#include "system_includes.h"
#include "add_suffix.h"
void  add_suffix(char *output_file, const char *suffix) {
  /*
    Catenate .suffix to output_file
    Called by: init_blue_combos, init_green_combos, init_red_combos
    Calls:     strlen, strcpy
    NB User should verify output file points to an buffer with enough
    space to add the suffix.
  */
  char *suff_pos;
  int out_len;
  out_len = strlen(output_file);
  output_file[out_len] = '.';
  output_file[out_len+1] = '\0';
  suff_pos = (char*)&output_file[out_len+1];
  strcpy(suff_pos,suffix);
}
