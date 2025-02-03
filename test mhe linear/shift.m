function [t0, x0, u0, con0] = shift(T, t0, x0, u,f, con_cov)
    st = x0;
    con0 = u(1,:)'+ sqrt(con_cov)*randn(2,1); %Applying noise to control input
    f_value = f(st,con0);
    st = st+ (T*f_value);
    x0 = full(st);

    t0 = t0 + T;
    u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end