    function [A,B,z0_block] = f_func_lin(z0, params,nStates)
M=params.magnet.m;
J=params.magnet.I(1);
mu0=params.physical.mu0;
Rad=params.solenoids.x(1);
rl=params.magnet.r; %Radius svevemagnet
rp=params.permanent.r(1); %Radius permanent magnet
ll=params.magnet.l; %Lengde svevemagnet
lp=params.permanent.l(1); %Lengde permanent magnet
ml = abs(params.magnet.J)*pi*rl^2*ll/mu0; % magnetisk moment svevemagnet
m = abs(params.permanent.J)*pi*rp^2*lp/mu0; % magnetisk moment permanent magnet
z_eq=z0;


q  = Rad^2 + z_eq^2;
p1 = -(4*Rad^4 - 27*Rad^2*z_eq^2 + 4*z_eq^4)/2;
p2 = -(z_eq*(11*Rad^4 - 23*Rad^2*z_eq^2 + z_eq^4))/4;
p3 = -10*(Rad*(Rad^2 - 4*z_eq^2))/2;
p4 = 2*(-3*Rad^4 + 24*Rad^2*z_eq^2 - 8*z_eq^4);
p5 = 1e3*(z_eq*(3*Rad^2 - 2*z_eq^2))/2;
p6 = 20*(z_eq*(4*Rad^2 - z_eq^2));
p7 = 20*(-Rad^4 + 13*Rad^2*z_eq^2 - z_eq^4);
p8 = 200*(Rad*z_eq)/2;

A =               [0,0,0,                                      1,  0,   0;
                   0,0,0,                                      0,  1,   0;
                   0,0,0,                                      0,  0,   1;
                   mu0*ml/(2*pi)*3*ml/M * p1/q^(9/2),0,-mu0*ml/(2*pi)*3*ml/M * p2/q^(9/2), 0 , 0 ,  0;
                   0 , -mu0*ml/(2*pi)*3*ml/M * p4/q^(9/2), 0 ,               0 , 0 ,  0;
                   -mu0*ml/(2*pi)*3*ml/J * p6/q^(7/2), 0 ,mu0*ml/(2*pi)* m/J * p7/q^(7/2), 0 , 0 ,  0];

B = mu0*ml/(2*pi)*[0,0;
                   0,0;
                   0,0;
                   3*mu0/M * p3/q^(7/2), 0;
                   0, 3*mu0/M * p5/q^(7/2)
                   -3*mu0/M * p8/q^(7/2), 0];

z0_block = zeros(nStates, 1);
z0_block(5) = mu0 * ml / (2 * pi) * (3 * ml / M * p4 / q^(9/2)) * z0;  % 5th row corresponds to z-dynamics

end