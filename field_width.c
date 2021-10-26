#include "system_includes.h"
#include "field_width.h"
int field_width(char *field, int length) {
  /*
    Determine the number of characters before the first , \n or \r in
    field, a string of length length
    Called by: csv_read
  */
  int iw;
  int not_sep;
  int i;
  int c;
  int comma;
  int lf;
  int cr;
  iw = 0;
  not_sep = 1;
  comma = ',';
  lf    = '\n';
  cr    = '\r';
  for (i=0;((i<length) && not_sep);i++) {
    c = field[i];
    not_sep = ((c != comma) && (c != lf) && (c!= cr));
    if (not_sep) iw += 1;
  }
  return(iw);
}
