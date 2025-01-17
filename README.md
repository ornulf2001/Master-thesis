For test mhe+mpc example:
1. Configure mpc_setup -> mhe_setup (model, weights, constraints, initial values, setpoint values)
2. Run either mpc alone (run_mpc.m) or mpc+mhe (run_mhe.m)

shift.m is responsible for applying control input to the model and moving to the next step.
plotting.m is for plotting states and inputs (nominal vs est.)
Draw_MPC_PS_Obstacles.m is for generating movie/GIF of the mpc in action.
