%+-----------------------------------------------------------------------+
%| Description:                                                          |
%| To derive the stiffness matrix in 2D space                            |
%+-----------------------------------------------------------------------+
%+-----------------------------------------------------------------------+
%| Created by Professor Hai Qing                                         |
%| Copyright: Nanjing University of Aeronautics and Astronautics 2016    |
%| Contact: qinghai@nuaa.edu.cn                                          |
%+-----------------------------------------------------------------------+
%| variables                                                             |
%| SF: shape function                                                    |
%| u: displacement                                                       |
%+-----------------------------------------------------------------------+
% function [K]=beam2
clc
syms L x EI EA c0 c1 c2 c3 u1 u2 v1 dv1 v2 dv2 R cosp sinp u1G v1G u2G v2G%real
EA=R*EI;
um=[v1 dv1 v2 dv2];
%% beam element
% eta=x/a;
v=c0+c1*x+c2*x^2+c3*x^3;  % deflection
dv=diff(v,x);
% boundary conditions
eq1=subs(v, x, 0)==v1;       
eq2=subs(dv, x, 0)==dv1;       
eq3=subs(v, x, L)==v2;
eq4=subs(dv, x, L)==dv2;       
% solve the unknown coefficients
s=solve([eq1,eq2,eq3,eq4],[c0,c1,c2,c3]);
% substitute the coefficient back into deflection function
% a=s.c0;
% b=s.c1;
% c=s.c2;
% d=s.c3;

un=s.c0+s.c1*x+s.c2*x^2+s.c3*x^3;
% reagrange and get the Shape Function: SF
SF=[s.c0,s.c1,s.c2,s.c3];
for kp=1:length(um)
    C=coeffs(un,um(kp));
    SF(kp)=C(2);   
end

v1L=-u1G*sinp+v1G*cosp;
v2L=-u2G*sinp+v2G*cosp;

Ubeam=EI/2*int((diff(SF(1),x,2)*v1L+diff(SF(2),x,2)*dv1+diff(SF(3),x,2)*v2L+diff(SF(4),x,2)*dv2)^2,x,0,L);

%% truss element
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
Ubar=L*EA/2*(diff(f1,x)*u1L+diff(f2,x)*u2L)^2;

%% frame element
U=Ubeam+Ubar;

k11=diff(U,u1G,u1G);
k12=diff(U,u1G,v1G);
k13=diff(U,u1G,dv1);
k14=diff(U,u1G,u2G);
k15=diff(U,u1G,v2G);
k16=diff(U,u1G,dv2);

k21=diff(U,v1G,u1G);
k22=diff(U,v1G,v1G);
k23=diff(U,v1G,dv1);
k24=diff(U,v1G,u2G);
k25=diff(U,v1G,v2G);
k26=diff(U,v1G,dv2);

k31=diff(U,dv1,u1G);
k32=diff(U,dv1,v1G);
k33=diff(U,dv1,dv1);
k34=diff(U,dv1,u2G);
k35=diff(U,dv1,v2G);
k36=diff(U,dv1,dv2);

k41=diff(U,u2G,u1G);
k42=diff(U,u2G,v1G);
k43=diff(U,u2G,dv1);
k44=diff(U,u2G,u2G);
k45=diff(U,u2G,v2G);
k46=diff(U,u2G,dv2);

k51=diff(U,v2G,u1G);
k52=diff(U,v2G,v1G);
k53=diff(U,v2G,dv1);
k54=diff(U,v2G,u2G);
k55=diff(U,v2G,v2G);
k56=diff(U,v2G,dv2);

k61=diff(U,dv2,u1G);
k62=diff(U,dv2,v1G);
k63=diff(U,dv2,dv1);
k64=diff(U,dv2,u2G);
k65=diff(U,dv2,v2G);
k66=diff(U,dv2,dv2);

stiffMatrix=simplify([k11 k12 k13 k14 k15 k16;k21 k22 k23 k24 k25 k26;k31 k32 k33 k34 k35 k36;...
    k41 k42 k43 k44 k45 k46;k51 k52 k53 k54 k55 k56;k61 k62 k63 k64 k65 k66]*(L/EI))