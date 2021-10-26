#ifndef _PROCESS_COMBO_H_
#define _PROCESS_COMBO_H_ 1
extern void process_combo(struct state_struct *state,
			  struct loop_state_struct *loop_state, 
			  int nzc,
			  int64_t *num_in_tca_p,
			  int *istep_index,
			  double x,
			  double y,
			  double z,
			  double sqdist);
#endif
