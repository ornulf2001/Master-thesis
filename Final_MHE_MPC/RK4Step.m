function x_next = RK4Step(f,x,U,dt,params)
    k1 = f(x,U,params);
    k2 = f(x+0.5*dt*k1,U,params);
    k3 = f(x+0.5*dt*k2,U,params);
    k4 = f(x + dt*k3, U, params);

    x_next = x+(dt/6)*(k1+2*k2+2*k3+k4);
end