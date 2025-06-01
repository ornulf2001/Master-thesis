Newest and most relevant release: **Final_MHE_MPC**
  - main.m runs the online control loop
  - MHE and MPC classes are defined in MHEClass.m and MPCclass.m, respectively. 




OLD RELEASES: **test mhe MSD**:
  - Here I develop an MHE routine for estimating the position and velocity of a Mass-spring-damper (MSD).
  - It should be robust for any other linear system, supply only A, B and C, and alter the dimensions/values of the weight matrices R, Q and M.
  - See /notes mhe setup.pdf for explanation of the cost function and constraints setup + dimensions
  - In /msd_sim.m we set up the simulator. Define initial value and control input sequence.
  
 
 
 

For **test mhe+mpc**: *Based on Mohamed W. Mehrez - [MPC and MHE implementation in Matlab using Casadi](https://www.youtube.com/watch?v=RrnkPrcpyEA&list=PLK8squHT_Uzej3UCUHjtOtm5X7pMFSgAL)*
1. Configure mpc_setup -> mhe_setup (model, weights, constraints, initial values, setpoint values)
2. Run either mpc alone (run_mpc.m) or mpc+mhe (run_mhe.m)

- shift.m is responsible for applying control input to the model and moving to the next step.
- plotting.m is for plotting states and inputs (nominal vs est.)
- Draw_MPC_PS_Obstacles.m is for generating movie/GIF of the mpc in action.
