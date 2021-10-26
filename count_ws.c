#include "system_includes.h"
#include "count_ws.h"
int count_ws(char *line) {
  /*
    Return the number of leading white space characters (<=32) 
    int the character string, line.\
    Called by: read_params
    Calls:     strlen(intrinsic).
  */
  int ws_chars;
  int i;
  int c;
  int line_len;
  /*
    Want to return 0 if line has zero length.
  */
  ws_chars = 0;
  line_len = strlen(line);
  for (i=0;i<line_len;i++) {
    c  = (int)line[i];
    if (c > 32) break;
    ws_chars += 1;
  }
  return (ws_chars);
}
