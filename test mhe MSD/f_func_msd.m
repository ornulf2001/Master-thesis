function [A,B,C,D]=f_func_msd(m,d,k)
A=[0,       1;
   -k/m, -d/m];
B=[0;1/m];
C=[1,0];
D=0;
end