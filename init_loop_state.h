#ifndef _INIT_LOOP_STATE_H_
#define _INIT_LOOP_STATE_H_ 1
extern int init_loop_state(struct state_struct *state,
			   int num_channels, 
			   double u_target, 
			   double v_target,
			   double test_rad);
#endif
