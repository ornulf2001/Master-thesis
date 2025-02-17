%%%%%%% Plotting %%%%%%%%%

figure(1)
%plot(Y_noisy(1,N_MHE+1:end)); hold on
plot(xsol2(1,:),"r"); hold on
plot(Y_sim(1,N_MHE+1:end),"k")

title("Estimated position. N: "+num2str(N_MHE)+", R: "+num2str(round(R,3))+", Q: "+mat2str(round(Q,3))+", M: "+mat2str(round(M,3)))
grid on
%legend("meas")
legend(["est. class","sim"])
xlabel("Time")
ylabel("x1")

figure(2)
plot(xsol2(2,:)); hold on
plot(X_sim(2,N_MHE+1:end))

legend(["est. class","sim"])
title("Estimated velcocity")
xlabel("Time")
ylabel("x2")