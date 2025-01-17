function plotting(N_MHE, mheIter, t,t_stop, xx, u_cl, st_est, u_est, con)
disp("Note that at linear speed -> 0, the system becomes unobservable and thus omega can't be estimated")


%%% Control inputs
figure(1)
subplot(211); hold on; grid on
ylabel('v (m/s)')
stairs(t(2:t_stop),u_cl(1:t_stop-1,1),'k','linewidth',1.5); axis([0 max(t) -2 2])
stairs(t(2:t_stop),con(1:t_stop-1,1),'g','linewidth',1); axis([0 max(t) -2 2])
stairs(t(N_MHE:N_MHE+mheIter-1),u_est(1:mheIter,1),'r','linewidth',1.5); axis([0 max(t) -2 2])


subplot(212); hold on; grid on
xlabel('time (s)')
ylabel('\omega (rad/s)')
stairs(t(2:t_stop),u_cl(1:t_stop-1,2),'k','linewidth',1.5); axis([0 max(t) -pi pi])
stairs(t(2:t_stop),con(1:t_stop-1,2),'g','linewidth',1); axis([0 max(t) -pi pi]) 
stairs(t(N_MHE:N_MHE+mheIter-1),u_est(1:mheIter,2),'r','linewidth',1.5); axis([0 max(t) -pi pi])


%%% States
figure(2)
subplot(311); hold on; grid on
ylabel('x [m]')
plot(t(1:t_stop),xx(1,1:t_stop),'k','linewidth',1.5); axis([0 max(t) -3 3]) 
plot(t(N_MHE:N_MHE+mheIter-1),st_est(1:mheIter,1),'r','linewidth',1.5); axis([0 max(t) -3 3])

subplot(312); hold on; grid on
ylabel('y [m]')
plot(t(1:t_stop),xx(2,1:t_stop),'k','linewidth',1.5); axis([0 max(t) -3 3])
plot(t(N_MHE:N_MHE+mheIter-1),st_est(1:mheIter,2),'r','linewidth',1.5); axis([0 max(t) -3 3])

subplot(313); hold on; grid on
ylabel('\theta [rad]')
xlabel('time (s)')
plot(t(1:t_stop),xx(3,1:t_stop),'k','linewidth',1.5); axis([0 max(t) -5 5])
plot(t(N_MHE:N_MHE+mheIter-1),st_est(1:mheIter,3),'r','linewidth',1.5); axis([0 max(t) -5 5])
