#include "system_includes.h"
#include "count_nws.h"
int count_nws(char *line) {
  /*
    Return the number of leading non-white space characters (<=32) 
    in the character string, line.\
    Called by: read_params
    Calls      strlen(intrinsic).
  */
  int nws_chars;
  int i;
  int c;
  int line_len;
  /*
    Want to return 0 if line has zero length.
  */
  nws_chars = 0;
  line_len = strlen(line);
  for (i=0;i<line_len;i++) {
    c  = (int)line[i];
    if (c <= 32) break;
    nws_chars += 1;
  }
  return (nws_chars);
}
