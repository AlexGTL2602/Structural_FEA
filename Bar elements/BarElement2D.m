clear all
close all
warning off
clc
syms x L EA c0 c1 cosp sinp u1 u2 u1G v1G u2G v2G 
%
% node infor (x1 y1)  (x2 y2)
% L=norm((x1-x2 y1-y2)) 
% cosp=(x1-x2)/L  sinp=(y1-y2)/L
%
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

u1L=u1G*cosp+v1G*sinp;
u2L=u2G*cosp+v2G*sinp;
U=L*EA/2*(diff(f1,x)*u1L+diff(f2,x)*u2L)^2;

k11=diff(U,u1G,u1G);
k12=diff(U,u1G,v1G);
k13=diff(U,u1G,u2G);
k14=diff(U,u1G,v2G);

k21=diff(U,v1G,u1G);
k22=diff(U,v1G,v1G);
k23=diff(U,v1G,u2G);
k24=diff(U,v1G,v2G);

k31=diff(U,u2G,u1G);
k32=diff(U,u2G,v1G);
k33=diff(U,u2G,u2G);
k34=diff(U,u2G,v2G);

k41=diff(U,v2G,u1G);
k42=diff(U,v2G,v1G);
k43=diff(U,v2G,u2G);
k44=diff(U,v2G,v2G);

stiffMatrix=[k11 k12 k13 k14;k21 k22 k23 k24;k31 k32 k33 k34;k41 k42 k43 k44]