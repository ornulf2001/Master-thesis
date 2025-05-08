%Bounds MPC
lbx = -inf;
ubx = inf;
lby = -inf;
uby = inf;
lbz = -inf;
ubz = inf;
lbphi = -inf;
ubphi = inf;
lbtheta = -inf;
ubtheta = inf;

lbxdot = -inf;
ubxdot = inf;
lbydot = -inf;
ubydot = inf;
lbzdot = -inf;
ubzdot = inf;
lbphidot = -inf;
ubphidot = inf;
lbthetadot = -inf;
ubthetadot = inf;

lbux = -2;
ubux = 2;
lbuz = -2;
ubuz = 2;

lbX1 = [lbx; lby; lbz; lbphi; lbtheta; lbxdot; lbydot; lbzdot; lbphidot; lbthetadot];
ubX1 = [ubx; uby; ubz; ubphi; ubtheta; ubxdot; ubydot; ubzdot; ubphidot; ubthetadot];
lbX = repmat(lbX1, N_MPC+1, 1);
ubX = repmat(ubX1, N_MPC+1, 1);

lbU1 = [lbux; lbuz];
ubU1 = [ubux; ubuz];
lbU = repmat(lbU1, N_MPC, 1);
ubU = repmat(ubU1, N_MPC, 1);

lb = [lbX; lbU];
ub = [ubX; ubU];

% Bounds uRef

lbuRef = [-inf; -inf; -inf; -inf;];
ubuRef = [inf; inf; inf; inf;];