#include "system_includes.h"
#include "loop_state_struct.h"
#include "gen_y_combo.h"
int gen_y_combo(struct loop_state_struct *loop_state, double *ycombo) {

  /*
    Get the next viable combination of y (green) channels.
    Called by: 
    Returns 1 if a new blue channel combination was found, 0 otherwise.


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
				    

    Fields of loop_state that are set:

       istep_index,         array of index values (length is num_channels)
                            On output the index values for the next
		            z(blue) channel combination.

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

  double delta_x;
  double delta_y;
  double delta_z;
  double z_new;
  double *vb_plane;
  double *ub_plane;

  int *istep_index;
  int *msteps_per_channel;
  int *is_xyz_p;

  int y_combos;
  int iy;

  int i3y;
  int i3step;

  int iiy;
  int nz;

  int istep;
  int ny;

  int inc1;
  int ithree;

  inc1                = 1;
  ithree              = 3;
  nz                  = loop_state->nz;
  ny                  = loop_state->ny;
  stepped_channel_xyz = loop_state->stepped_channel_xyz;
  prev_xyzs           = loop_state->prev_xyz;
  /*
  curr_xyz            = loop_state->curr_xyz;
  */
  istep_index         = loop_state->istep_index;
  msteps_per_channel  = loop_state->msteps_per_channel;
  is_xyz_p            = loop_state->istepped_xyz_ptrs;
  vb_plane            = loop_state->vb_plane;
  ub_plane            = loop_state->ub_plane;

  y_combos = 0;
  for (iy=nz+ny-1;((iy>=nz) && (y_combos == 0));iy--) {
    i3y = iy + iy + iy;
    if (istep_index[iy] < msteps_per_channel[iy]-1) {
      y_combos = 1;
      istep_index[iy] += 1;
      istep = istep_index[iy];
      i3step = istep + istep + istep;
      channel_steps = &stepped_channel_xyz[is_xyz_p[iy]];
      prev_xyz = &prev_xyzs[i3y];
      delta_z  = channel_steps[i3step+2];
      delta_y  = channel_steps[i3step+1];
      delta_x  = channel_steps[i3step];
      /*
      curr_xyz[0] = prev_xyz[0] + delta_x;
      curr_xyz[1] = prev_xyz[1] + delta_y;
      curr_xyz[2] = prev_xyz[2] + delta_z;
      */
      ycombo[0] = prev_xyz[0] + delta_x;
      ycombo[1] = prev_xyz[1] + delta_y;
      ycombo[2] = prev_xyz[2] + delta_z;
      /*
      ycombo[3] = ddot_(&ithree,ycombo,&inc1,vb_plane,&inc1);
      ycombo[4] = ddot_(&ithree,ycombo,&inc1,ub_plane,&inc1);
      */
      ycombo[3] = ((ycombo[0] * vb_plane[0]) + (ycombo[1] * vb_plane[1])) +
	           (ycombo[2] * vb_plane[2]);
      ycombo[4] = ((ycombo[0] * ub_plane[0]) + (ycombo[1] * ub_plane[1])) +
	           (ycombo[2] * ub_plane[2]);
      /*
	Propagate forward curr_xyz into prev_xyz for later
	y channels
      */
      prev_xyzf = prev_xyz + 3; /* Caution address arithmetic */
      for (iiy = iy+1;iiy < nz+ny;iiy++) {
	/*
	prev_xyzf[0] = curr_xyz[0];
	prev_xyzf[1] = curr_xyz[1];
	prev_xyzf[2] = curr_xyz[2];
	*/
	prev_xyzf[0] = ycombo[0];
	prev_xyzf[1] = ycombo[1];
	prev_xyzf[2] = ycombo[2];
	prev_xyzf    = prev_xyzf + 3; /* Caution address arithmetic */
      }
    } else {
      istep_index[iy] = 0;
    } /* end else channel iy is in range.*/
  } /* end for (iy ...) */
  return(y_combos);
}
