%%%%%%%
%This file is used to provide the parameters about vehicle 
%%%%%%%
m = 1573;
Iz = 2873;
l_f = 1.35;
l_r = 1.35;
W_f=l_r/(l_f+l_r);
W_r=l_f/(l_f+l_r);
p=[m, Iz, l_f,l_r];
Cm1=17303;
Cm2=175;
Cr0=120;
Cr2=0.5*1.225*0.35*2.5;
B_r = 13;
C_r = 2;
D_r = W_f*m*9.81*1.2;
B_f = 13;
C_f = 2;
D_f = W_r*m*9.81*1.2;
pf=[Cm1,Cm2,Cr0,Cr2,B_r,C_r,D_r,B_f,C_f,D_f];
