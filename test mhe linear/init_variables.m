import casadi.*

T = 0.01;
z_eq=0.0365;
x_eq = [0,z_eq,0,0,0,0]';
u_eq = [0,0]';

x = SX.sym('x'); 
z = SX.sym('z'); 
theta = SX.sym('theta');
xdot = SX.sym('xdot'); 
zdot = SX.sym('zdot'); 
thetadot = SX.sym('thetadot');
states = [x; z; theta; xdot; zdot; thetadot]; 
nStates = length(states);

ux = SX.sym('ux'); 
uz = SX.sym('uz');
controls = [ux;uz]; 
nControls = length(controls);

params=getParams();

function data = getParams()
    persistent params

    if isempty(params)
        parameters;
    end
    data = params;
end
