#include "system_includes.h"
#include "extract_root.h"
void extract_root(char *input_file, char *output_file) {
  /*
    Remove the .* suffix from input_file putting the null_terminated
    resulting string in output_file
    Called by: init_blue_combos, init_green_combos, init_red_combos
    Calls:     strlen,strncpy
  */
  int root_len;
  int in_len;
  int i;
  int ipad;
  in_len = strlen(input_file);
  root_len = in_len;
  for (i=in_len-1;i>= 0;i--) {
    if (input_file[i] == '.') {
      root_len = i;
      break;
    }
  }
  strncpy(output_file,input_file,root_len);
  output_file[root_len] = '\0';
}
