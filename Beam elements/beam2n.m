
%+-----------------------------------------------------------------------+
%| variables                                                             |
%| SF: shape function                                                    |
%| u: displacement                                                       |
%+-----------------------------------------------------------------------+
function [K]=beam2n
clc
syms a b c d eta L x EI c0 c1 c2 c3 v1 dv1 v2 dv2 cosA sinA %real
um=[v1 dv1 v2 dv2];
% eta=x/(L/2);
v=c0+c1*eta+c2*eta^2+c3*eta^3;  % deflection
dv=diff(v,eta)/(L/2);
% boundary conditions
eq1=subs(v, eta, -1)==v1;       
eq2=subs(dv, eta, -1)==dv1;       
eq3=subs(v, eta, 1)==v2;
eq4=subs(dv, eta, 1)==dv2;       
% solve the unknown coefficients
s=solve([eq1,eq2,eq3,eq4],[c0,c1,c2,c3]);
% substitute the coefficient back into deflection function
% a=s.c0;
% b=s.c1;
% c=s.c2;
% d=s.c3;

un=s.c0+s.c1*eta+s.c2*eta^2+s.c3*eta^3;
% reagrange and get the Shape Function: SF
SF=[s.c0,s.c1,s.c2,s.c3];
for kp=1:length(um)
    C=coeffs(un,um(kp));
    SF(kp)=C(2);   
end
SF2=diff(SF,eta,2)/(L/2)^2;
% calculate the stiffness matrix
Ke=int(EI*transpose(SF2)*SF2,eta,-1,1)*(L/2)
