%%%%%%% Plotting %%%%%%%%%
figure(1)
%plot(Y_noisy(1,N_MHE+1:end),"k") 
hold on
%plot(xsol(1,:),"r"); hold on
%plot(Y_noisy(1,N_MHE+1:end))
plot(xsol2(1,:),"r");
plot(Y_sim(1,N_MHE+1:end),"k")
title("Estimated x")%. N: "+num2str(N_MHE)+", R: "+num2str(round(R,3))+", Q: "+mat2str(round(Q,3))+", M: "+mat2str(round(M,3)))
grid on
%legend("meas")
legend(["est. class","sim"])
xlabel("Time")
ylabel("x")
%ylim([-1,1]);

figure(2)
%plot(xsol(2,:)); hold on
plot(xsol2(2,:));
hold on
plot(Y_sim(2,N_MHE+1:end))
legend(["est. class","sim"])
title("Estimated z")
xlabel("Time")
ylabel("z")

figure(3)
plot(xsol2(3,:));
hold on
plot(X_sim(3,N_MHE+1:end))
legend(["est. class","sim"])
title("Estimated Theta")
xlabel("Time")
ylabel("theta")

figure(4)
hold on; grid on;
plot(xsol2(1,:),xsol2(2,:))
xlabel("x")
ylabel("z")