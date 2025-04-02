function y = maglevSystemMeasurements2d(x,u,params)
    % ### HACK from 2d to 3d ###
    x_full = [
        x(1),0,x(2),0,x(3),0,...
        x(4),0,x(5),0,x(6),0
        ]';
    u_full = [u(1), u(2), 0, 0]';
    % ##########################
    
    y = maglevSystemMeasurements(x_full,u_full,params);
    
    % ### HACK from 3d to 2d ###
    y = y([1,3]');
    % ##########################
end