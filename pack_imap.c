#include "system_includes.h"
#include "pack_imap.h"
void pack_imap(int num_channels,int *imap, int *imap_packed) {
  /*
    Pack the imap array into bytes of imap_packed packing into a
    eight bit fields.
    Called by: init_z_dom_combos, init_y_dom_combos, init_x_dom_combos
  */
  int left_shift;
  int packed_entry;

  int packed_pos;
  int mask;

  int new_piece;
  int i;

  left_shift = 24;
  packed_entry = 0;
  packed_pos   = 0;
  mask         = 255;
  for (i=0;i<num_channels;i++) {
    new_piece = imap[i] & mask;
    if (left_shift > 0) {
      new_piece = new_piece << left_shift;
      packed_entry = packed_entry | new_piece;
      left_shift = left_shift - 8;
    } else {
      packed_entry = packed_entry | new_piece;
      imap_packed[packed_pos] = packed_entry;
      left_shift = 24;
      packed_entry = 0;
      packed_pos = packed_pos + 1;
    }
  }
  if (left_shift != 24) {
    imap_packed[packed_pos] = packed_entry;
  }
}
