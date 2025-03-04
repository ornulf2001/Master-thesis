function rhs=f_func(states,controls,params)
    import casadi.*
    x=states(1); z=states(2); theta=states(3);
    xdot=states(4); zdot=states(5); thetadot=states(6);
    ux=controls(1); uz=controls(2);

    %params 
    mu0=params.physical.mu0;
    gr=params.physical.g;
    
    Rad=params.solenoids.x(1); %Avstand fra nullpunkt til solenoid
    rl=params.magnet.r; %Radius svevemagnet
    rp=params.permanent.r(1); %Radius permanent magnet
    ll=params.magnet.l; %Lengde svevemagnet
    lp=params.permanent.l(1); %Lengde permanent magnet
    
    mass=params.magnet.m; %Mass
    J=params.magnet.I(1); %Moment of inertia
    ml = abs(params.magnet.J)*pi*rl^2*ll/mu0; % magnetisk moment svevemagnet
    m = abs(params.permanent.J)*pi*rp^2*lp/mu0; % magnetisk moment permanent magnet
    ml_vec=[0;0;ml];
    
    % hybridmagneter
    Ip = (ux + uz)/2; %Strøm positiv hybrid magnet
    In = (uz - ux)/2; %Strøm negativ hybrid magnet
    mp= (m+mu0*Ip).*[sin(theta);0;-cos(theta)];
    mn= (m+mu0*In).*[sin(theta);0;-cos(theta)];
    rn = [x - Rad*cos(theta);
          0;
          z - Rad*sin(theta)];
    rp = [x + Rad*cos(theta);
          0;
          z + Rad*sin(theta)];
    
    %Magnetisk skalarpotensiale, Phi
    phi=mu0/(4*pi) * (dot(mp,rp)/norm(rp)^3 + dot(mn,rn)/norm(rn)^3);
    
    dphi_dx = jacobian(phi, x);
    d2phi_dxdz = jacobian(dphi_dx, z);
    dphi_dz = jacobian(phi, z);
    d2phi_dz2 = jacobian(dphi_dz, z);
    
    xddot = -(ml/mass) * d2phi_dxdz;
    zddot = -(m/mass) * d2phi_dz2 - gr;
    thetaddot = -(ml/J) * dphi_dx;
   

    rhs=[xdot;zdot;thetadot;xddot;zddot;thetaddot];
end