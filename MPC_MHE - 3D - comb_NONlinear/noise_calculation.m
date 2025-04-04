clc,clear

load("data_no_control 1.mat")
Y_noisy=[data.y.bx0';data.y.by0';-data.y.bz0';data.y.bx1';data.y.by1';-data.y.bz1';data.y.bx2';data.y.by2';-data.y.bz2']*1e-3;
R=inv(cov(Y_noisy'));
save("noise_cov","R")