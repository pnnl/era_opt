#include "system_includes.h"
#include "pack_isteps.h"
void  pack_isteps(int num_steps, int sbpi, int *isteps, int *isteps_packed) {
  /*
    Pack the imap array into bytes of imap_packed packing into a
    eight bit fields.
  */
  int start_left_shift;
  int left_shift;

  int packed_entry;
  int packed_pos;

  int mask;
  int new_piece;

  int i;
  int shift_delta;

  if (sbpi == 8) {
    start_left_shift = 24;
    mask       = 255;
    shift_delta = 8;
  } else {
    start_left_shift = 16;
    mask       = 65535;
    shift_delta = 16;
  }
  left_shift = start_left_shift;
  packed_entry = 0;
  packed_pos   = 0;
  for (i=0;i<num_steps;i++) {
    new_piece = isteps[i] & mask;
    if (left_shift > 0) {
      new_piece = new_piece << left_shift;
      packed_entry = packed_entry | new_piece;
      left_shift = left_shift - shift_delta;
    } else {
      packed_entry = packed_entry | new_piece;
      isteps_packed[packed_pos] = packed_entry;
      left_shift = start_left_shift;
      packed_entry = 0;
      packed_pos = packed_pos + 1;
    }
  }
  if (left_shift != start_left_shift) {
    isteps_packed[packed_pos] = packed_entry;
  }
}
