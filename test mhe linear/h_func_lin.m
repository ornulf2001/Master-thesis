function [C,D] = h_func_lin(z0,params)
M=params.magnet.m;
mu0=params.physical.mu0;

C=[3*M*mu0/(4*pi*z0^4), 0, M*mu0/(4*pi*z0^3),0,0,0;
   0,-3*M*mu0/(2*pi*z0^4),0,0,0,0];

D= mu0^2/(2*pi*z0^4)* [0,0;
                       1,1];

end