#include "system_includes.h"
#include "pos_char.h"
int pos_char(char *value, char ic) {
  /*
    Return the position (starting from 1 of a character, ic, in a string,
    0 if ic is not in the string.
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
    if (*string == ic) {
      pos = i+1;
      break;
    }
    string++; /* Caution address arithmetic */
  }
  return(pos);
}
