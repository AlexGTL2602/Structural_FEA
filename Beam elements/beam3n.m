
%+-----------------------------------------------------------------------+
%| variables                                                             |
%| SF: shape function                                                    |
%| u: displacement                                                       |
%+-----------------------------------------------------------------------+
function [K]=beam3n
clc
syms L EI c0 c1 c2 c3 c4 c5 eta v1 dv1 v2 dv2 v3 dv3 cosA sinA %real
um=[v1 dv1 v2 dv2 v3 dv3];
% eta=x/a;
v=c0+c1*eta+c2*eta^2+c3*eta^3+c4*eta^4+c5*eta^5;  % deflection
dv=diff(v,eta)/(L/2);
% boundary conditions
eq1=subs(v, eta, -1)==v1;       
eq2=subs(dv, eta, -1)==dv1;       
eq3=subs(v, eta, 0)==v2;
eq4=subs(dv, eta, 0)==dv2;    
eq5=subs(v, eta, 1)==v3;
eq6=subs(dv, eta, 1)==dv3; 
% solve the unknown coefficients
s=solve([eq1,eq2,eq3,eq4,eq5,eq6],[c0,c1,c2,c3,c4,c5]);
% substitute the coefficient back into deflection function
% a=s.c0;
% b=s.c1;
% c=s.c2;
% d=s.c3;

un=s.c0+s.c1*eta+s.c2*eta^2+s.c3*eta^3+s.c4*eta^4+s.c5*eta^5;
% reagrange and get the Shape Function: SF
SF=[s.c0,s.c1,s.c2,s.c3,s.c4,s.c5];
for kp=1:length(um)
    C=coeffs(un,um(kp));
    SF(kp)=C(2);   
end
SF2=diff(SF,eta,2)/(L/2)^2;
% calculate the stiffness matrix
Ke=int (EI*transpose(SF2)*SF2,eta,-1,1)*L/2
