#include "system_includes.h"
#include "loop_state_struct.h"
#include "bf_gen_z_combo.h"
int bf_gen_z_combo(struct loop_state_struct *loop_state, double *zcombo) {

  /*
    Get the next viable combination of z dominant (blue) channels.
    Called by: gen_z_dom_combos
    Returns 1 if a new z dominant channel combination was found, 0 otherwise.


    Arguments :             TMF     Description
       loop_state           p*b     loop_state structure.

    Fields of loop_state that are used:
       nz                   number of z-dominant channels

       istep_index,         array of index values (length is num_channels)
                            On input the index values for the previous z combo.

       msteps_per_channel,  maximum number of steps per channel
                            For multiscale note that msteps_per_channel
			    has been set to max_index+1 in all channels.

       stepped_channel_xyz, 3 by total steps array with the scaled
                            x,y,z tristimulus values at each power level,
			    for each channel. Indexd by is_xyz_p.

       istepped_xyz_ptrs    Index vector of length nun_channels pointing to
                            the first element of the tristimulus values for
			    each channel, referred to as is_xyz_p internally

       prev_xyzs            (3 x num_channels),
                            on input the previous tristimulus values for
                            channels at the power levels specified in istep_index
				    
       zmax                 maximum z contribution.

    Fields of loop_state that are set:

       istep_index,         array of index values (length is num_channels)
                            On output the index values for the next
		            z dominant (blue) channel combination.

       prev_xyzs            (3 x num_channels),
                            on output the current  tristimulus values for
			    channels at the power levels specified in istep_index
				    

       
       curr_xyz             (x,y,z) values for new combination of blue channels

  */
  double *stepped_channel_xyz; 
  double *prev_xyzs;
  /*
  double *curr_xyz;
  */

  double *channel_steps;
  double *prev_xyz;
  double *prev_xyzf;

  double zmax; 

  double delta_x;
  double delta_y;
  double delta_z;
  double z_new;

  int *istep_index;
  int *msteps_per_channel;
  int *is_xyz_p;

  int z_combos;
  int iz;

  int i3z;
  int i3step;

  int iiz;
  int nz;

  int istep;
  int ipad;

  nz                  = loop_state->nz;
  stepped_channel_xyz = loop_state->stepped_channel_xyz;
  prev_xyzs           = loop_state->prev_xyz;
  /*
  curr_xyz            = loop_state->curr_xyz;
  */
  istep_index         = loop_state->istep_index;
  msteps_per_channel  = loop_state->msteps_per_channel;
  is_xyz_p            = loop_state->istepped_xyz_ptrs;
  zmax                = loop_state->zmax;

  z_combos = 0;
  for (iz=nz-1;((iz>=0) && (z_combos == 0));iz--) {
    i3z = iz + iz + iz;
    if (istep_index[iz] < msteps_per_channel[iz]-1) {
      istep_index[iz] += 1;
      istep = istep_index[iz];
      i3step = istep + istep + istep;
      channel_steps = &stepped_channel_xyz[is_xyz_p[iz]];
      prev_xyz = &prev_xyzs[i3z];
      delta_z  = channel_steps[i3step+2];
      z_new    = prev_xyz[2] + delta_z;
      z_combos = 1;
      delta_y  = channel_steps[i3step+1];
      delta_x  = channel_steps[i3step];
      /*
      curr_xyz[0] = prev_xyz[0] + delta_x;
      curr_xyz[1] = prev_xyz[1] + delta_y;
      curr_xyz[2] = z_new;
      */
      zcombo[0] = prev_xyz[0] + delta_x;
      zcombo[1] = prev_xyz[1] + delta_y;
      zcombo[2] = z_new;

      /*
        Propagate forward curr_xyz into prev_xyz for later
        z channels
      */
      prev_xyzf = prev_xyz + 3; /* Caution address arithmetic */
      for (iiz = iz+1;iiz < nz;iiz++) {
        /*
        prev_xyzf[0] = curr_xyz[0];
        prev_xyzf[1] = curr_xyz[1];
        prev_xyzf[2] = curr_xyz[2];
        */
        prev_xyzf[0] = zcombo[0];
        prev_xyzf[1] = zcombo[1];
        prev_xyzf[2] = zcombo[2];
        prev_xyzf    = prev_xyzf + 3; /* Caution address arithmetic */
      }
    } else {
      istep_index[iz] = 0;
    } /* end else iz channel has run all steps bump to earlier channel. */
  } /* end for (iz ...) */
  return(z_combos);
}
