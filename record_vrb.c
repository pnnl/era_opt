#include "system_includes.h"
#include "state_struct.h"
#include "pos_char.h"
#include "record_vrb.h"
void record_vrb(char *value, struct state_struct *state) {
  /*
    Read a commas separated list of first_index,last_index,radius
        triples that are preceded by the number of triples followed by a 
	colon.
    Called by: read_params
    Calls:     pos_char, sscanf,
  */
  double *vrb_radii;
  int    *vrb_first_index;
  int    *vrb_last_index;
  char *string;
  double radius;


  int number_vrb;
  int ns;

  int i;
  int inum;

  int first_index;
  int last_index;

  int loc_comma;
  int loc_colon;

  FILE *lfp;

  string = value;
  vrb_radii       = state->vrb_radii;
  vrb_first_index = state->vrb_first_index;
  vrb_last_index  = state->vrb_last_index;

  ns = sscanf(string,"%d",&number_vrb);
  if (ns == 1) {
    if (number_vrb <= 0) {
      state->number_vrb = 0;
      state->use_vrb    = 0;
    } else {
      loc_colon = pos_char(string,':');
      string += loc_colon; /* Caution address arithmetic. */
      inum = 0;
      for (i=0;i<number_vrb;i++) {
	ns = sscanf(string,"%d",&first_index);
	if (ns == 1) {
	  loc_comma = pos_char(string,',');
	  string += loc_comma; /* Caution address arithmetic. */
	  ns = sscanf(string,"%d",&last_index);
	  if (ns == 1) {
	    loc_comma = pos_char(string,',');
	    string += loc_comma; /* Caution address arithmetic. */
	    ns == sscanf(string,"%le",&radius);
	    if (ns == 1) {
	      loc_comma = pos_char(string,',');
	      string += loc_comma; /* Caution address arithmetic. */
	      vrb_radii[inum] = radius;
	      vrb_first_index[inum] = first_index;
	      vrb_last_index[inum]  = last_index;
	      inum += 1;
	    } else {
	      break;
	    }
	  } else {
	    break;
	  }
	} else {
	  break;
	}
      } /* end for (i...) */
      state->number_vrb = inum;
      state->use_vrb    = 1;
    }
  }
}
