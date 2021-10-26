#include "system_includes.h"
#include "upcase.h"
void upcase (int sl, char *string)  {
  /*
    Convert lower case to upper case in a string that may
    have non-alphabetic characters included.
    Called by: read_params
  */
  int i;
  int ic;
  for (i=0;i<sl;i++) {
    ic = (int)string[i];
    if ((ic > 96) && (ic < 123)) {
      ic = ic - 32;
    }
    string[i] = (char)ic;
  }
}

