clear all
close all
warning off
clc
syms x L EA c0 c1 u1 u2

u=c0+c1*x;
eq1=subs(u-u1,x,0);
eq2=subs(u-u2,x,L);
sol=solve([eq1==0,eq2==0],[c0 c1]);
C0=sol.c0;
C1=sol.c1;
ff=coeffs(subs(u,[c0 c1],[C0 C1]),u1);
f1=ff(2)
ff=coeffs(subs(u,[c0 c1],[C0 C1]),u2);
f2=ff(2)

U=L*EA/2*(diff(f1,x)*u1+diff(f2,x)*u2)^2;
k11=diff(U,u1,u1);
k12=diff(U,u1,u2);
k21=diff(U,u2,u1);
k22=diff(U,u2,u2);
stiffMatrix=[k11 k12;k21 k22]