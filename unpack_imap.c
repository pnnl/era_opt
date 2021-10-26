#include "system_includes.h"
#include "unpack_imap.h"
void unpack_imap(int num_channels,int *imap, int *imap_packed) {
  /*
    Unpack the imap_packed array into imap array with one element
    per map position.
    Called by:
  */
  int right_shift;
  int packed_entry;

  int packed_pos;
  int unpacked_pos;

  int mask;
  int new_entry;

  int i;
  int ipad;

  right_shift = 24;
  packed_entry = imap_packed[0];
  packed_pos   = 0;
  unpacked_pos = 0;
  mask         = 255;
  for (i=0;i<num_channels;i++) {
    if (right_shift > 0) {
      new_entry = (packed_entry >> right_shift) & mask;
      right_shift = right_shift - 8;
    } else {
      new_entry = packed_entry & mask;
      right_shift = 24;
      packed_pos += 1;
      if (i < num_channels - 1) {
	packed_entry = imap_packed[packed_pos];
      }
    }
    imap[unpacked_pos] = new_entry;
    unpacked_pos += 1;
  }
}
