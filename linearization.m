function [sys_d,Afun,Bfun]=linearization(state,input,pf,p,T)
%% following part is used to test this function 



%vecotr parameter of force
Cm1=pf(1);
Cm2=pf(2);
Cr0=pf(3);
Cr2=pf(4);
B_r=pf(5);
C_r=pf(6);
D_r=pf(7);
B_f=pf(8);
C_f=pf(9);
D_f=pf(10);
m=p(1);
Iz=p(2);
l_f=p(3);
l_r=p(4);
%
persistent A B x y phi v_x v_y omega delta d 
if isempty(A)
x=sym('x','real');
y=sym('y','real');
phi=sym('phi','real');
v_x=sym('v_x','real');
v_y=sym('v_y','real');
omega=sym('omega','real');
delta=sym('delta','real');
d=sym('d','real');
    alpha_f = -atan2(l_f*omega + v_y,v_x)+delta;
    alpha_r =  atan2(l_r*omega - v_y,v_x);

    F_fy = D_f*sin(C_f*atan(B_f*alpha_f));
    F_ry = D_r*sin(C_r*atan(B_r*alpha_r));

    F_rx = (Cm1*d-Cm2*d*v_x-Cr0-Cr2*v_x^2);

  f1=v_x*cos(phi) - v_y*sin(phi);
  f2=v_y*cos(phi) + v_x*sin(phi);
  f3=omega;
  f4=1/m*(F_rx - F_fy*sin(delta) + m*v_y*omega);
  f5=1/m*(F_ry + F_fy*cos(delta) - m*v_x*omega);
  f6=1/Iz*(F_fy*l_f*cos(delta)- F_ry*l_r);

 A=[diff(f1,x),diff(f1,y),diff(f1,phi),diff(f1,v_x),diff(f1,v_y),diff(f1,omega);
      diff(f2,x),diff(f2,y),diff(f2,phi),diff(f2,v_x),diff(f2,v_y),diff(f2,omega);
      diff(f3,x),diff(f3,y),diff(f3,phi),diff(f3,v_x),diff(f3,v_y),diff(f3,omega);
      diff(f4,x),diff(f4,y),diff(f4,phi),diff(f4,v_x),diff(f4,v_y),diff(f4,omega);
      diff(f5,x),diff(f5,y),diff(f5,phi),diff(f5,v_x),diff(f5,v_y),diff(f5,omega);
      diff(f6,x),diff(f6,y),diff(f6,phi),diff(f6,v_x),diff(f6,v_y),diff(f6,omega);
 ];
 B= [diff(f1,d),diff(f1,delta);
       diff(f2,d),diff(f2,delta);
       diff(f3,d),diff(f3,delta);
       diff(f4,d),diff(f4,delta);
       diff(f5,d),diff(f5,delta);
       diff(f6,d),diff(f6,delta) ];
end;

A_c= double(subs(A,{x,y,phi,v_x,v_y,omega,delta,d}, {state(1),state(2),state(3),state(4),state(5),state(6),input(1),input(2)} ));
B_c= double(subs(B,{x,y,phi,v_x,v_y,omega,delta,d}, {state(1),state(2),state(3),state(4),state(5),state(6),input(1),input(2)} ));
C_c=[1 0 0 0 0 0
    0 1 0  0 0 0
    0 0 1 0 0 0];
D_c=0;
sys_c= ss(A_c, B_c, C_c, D_c );
sys_d= c2d(sys_c,T);
Afun=A;
Bfun=B;


       
       
       
  
 
 
 
 
 
 