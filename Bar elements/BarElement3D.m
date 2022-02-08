clear all
close all
warning off
clc
syms x L EA c0 c1 l m n sinp u1 u2 u1G v1G w1G u2G v2G w2G
%
% node infor (x1 y1 z1)  (x2 y2 z2 )
% L=norm((x1-x2 y1-y2 z1-z2)) 
% l=(x1-x2)/L  m=(y1-y2)/L  n=(z1-z2)/L
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

u1L=u1G*l+v1G*m+w1G*n;
u2L=u2G*l+v2G*m+w2G*n;
U=L*EA/2*(diff(f1,x)*u1L+diff(f2,x)*u2L)^2;
k11=diff(U,u1G,u1G);
k12=diff(U,u1G,v1G);
k13=diff(U,u1G,w1G);
k14=diff(U,u1G,u2G);
k15=diff(U,u1G,v2G);
k16=diff(U,u1G,w2G);

k21=diff(U,v1G,u1G);
k22=diff(U,v1G,v1G);
k23=diff(U,v1G,w1G);
k24=diff(U,v1G,u2G);
k25=diff(U,v1G,v2G);
k26=diff(U,v1G,w2G);

k31=diff(U,w1G,u1G);
k32=diff(U,w1G,v1G);
k33=diff(U,w1G,w1G);
k34=diff(U,w1G,u2G);
k35=diff(U,w1G,v2G);
k36=diff(U,w1G,w2G);

k41=diff(U,u2G,u1G);
k42=diff(U,u2G,v1G);
k43=diff(U,u2G,w1G);
k44=diff(U,u2G,u2G);
k45=diff(U,u2G,v2G);
k46=diff(U,u2G,w2G);

k51=diff(U,v2G,u1G);
k52=diff(U,v2G,v1G);
k53=diff(U,v2G,w1G);
k54=diff(U,v2G,u2G);
k55=diff(U,v2G,v2G);
k56=diff(U,v2G,w2G);

k61=diff(U,w2G,u1G);
k62=diff(U,w2G,v1G);
k63=diff(U,w2G,w1G);
k64=diff(U,w2G,u2G);
k65=diff(U,w2G,v2G);
k66=diff(U,w2G,w2G);

stiffMatrix=[k11 k12 k13 k14 k15 k16;k21 k22 k23 k24 k25 k26;k31 k32 k33 k34 k35 k36;k41 k42 k43 k44 k45 k46;k51 k52 k53 k54 k55 k56;k61 k62 k63 k64 k65 k66]