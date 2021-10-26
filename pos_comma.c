#include "system_includes.h"
#include "pos_comma.h"
int pos_comma(char *value) {
  /*
    Return the position (starting from 1 of a comma in a string,
    0 if no comma in the string.
    Called by: record_steps_per_channel
    Calls:     strlen
  */
  char *string;
  int i;
  int pos;
  int len;
  int ipad;
  pos = 0;
  string = value;
  len    = strlen(value);
  for (i=0;i<len;i++) {
    if (*string == ',') {
      pos = i+1;
      break;
    }
    string++; /* Caution address arithmetic */
  }
  return(pos);
}
