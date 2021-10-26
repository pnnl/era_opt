#include "system_includes.h"
#include "unpack_isteps.h"
void  unpack_isteps(int num_steps, int sbpi, int *isteps, int *isteps_packed) {
  /*
    Unpack the imap_packed array into the  imap array.
  */
  int start_right_shift;
  int right_shift;

  int packed_entry;
  int packed_pos;

  int unpacked_pos;

  int mask;
  int new_entry;

  int i;
  int shift_delta;

  if (sbpi == 8) {
    start_right_shift = 24;
    mask       = 255;
    shift_delta = 8;
  } else {
    start_right_shift = 16;
    mask       = 65535;
    shift_delta = 16;
  }
  right_shift = start_right_shift;
  packed_entry = isteps_packed[0];
  packed_pos   = 0;
  unpacked_pos = 0;
  for (i=0;i<num_steps;i++) {
    if (right_shift > 0) {
      new_entry = (packed_entry >> right_shift) & mask;
      right_shift = right_shift - shift_delta;
    } else {
      new_entry = packed_entry & mask;
      right_shift = start_right_shift;
      packed_pos  += 1;
      if (i < (num_steps - 1)) {
	packed_entry = isteps_packed[packed_pos];
      }
    }
    isteps[unpacked_pos] = new_entry;
    unpacked_pos += 1;
  }
}
