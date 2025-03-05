%Bounds MPC
lbx = -inf;
ubx = inf;
lbz = -inf;
ubz = inf;
lbtheta = -inf;
ubtheta = inf;

lbxdot = -inf;
ubxdot = inf;
lbzdot = -inf;
ubzdot = inf;
lbthetadot = -inf;
ubthetadot = inf;

lbux = -inf;
ubux = inf;
lbuz = -inf;
ubuz = inf;

lbX1 = [lbx; lbz; lbtheta; lbxdot; lbzdot; lbthetadot];
ubX1 = [ubx; ubz; ubtheta; ubxdot; ubzdot; ubthetadot];
lbX = repmat(lbX1, N_MPC+1, 1);
ubX = repmat(ubX1, N_MPC+1, 1);

lbU1 = [lbux; lbuz];
ubU1 = [ubux; ubuz];
lbU = repmat(lbU1, N_MPC, 1);
ubU = repmat(ubU1, N_MPC, 1);

lb = [lbX; lbU];
ub = [ubX; ubU];
