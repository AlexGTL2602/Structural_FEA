%+-----------------------------------------------------------------------+
%| Description:                                                          |
%| To derive the stiffness matrix in 2D space                            |
%+-----------------------------------------------------------------------+
%| variables                                                             |
%| SF: shape function                                                    |
%| u: displacement                                                       |
%+-----------------------------------------------------------------------+
function [K]=beam3
clc
syms a b c d e f L x EI c0 c1 c2 c3 c4 c5 v1 dv1 v2 dv2 v3 dv3 cosA sinA %real
um=[v1 dv1 v2 dv2 v3 dv3];
% eta=x/a;
v=c0+c1*x+c2*x^2+c3*x^3+c4*x^4+c5*x^5;  % deflection
dv=diff(v,x);
% boundary conditions
eq1=subs(v, x, 0)==v1;       
eq2=subs(dv, x, 0)==dv1;       
eq3=subs(v, x, L/2)==v2;
eq4=subs(dv, x, L/2)==dv2;    
eq5=subs(v, x, L)==v3;
eq6=subs(dv, x, L)==dv3; 
% solve the unknown coefficients
s=solve([eq1,eq2,eq3,eq4,eq5,eq6],[c0,c1,c2,c3,c4,c5]);
% substitute the coefficient back into deflection function
% a=s.c0;
% b=s.c1;
% c=s.c2;
% d=s.c3;

un=s.c0+s.c1*x+s.c2*x^2+s.c3*x^3+s.c4*x^4+s.c5*x^5;
% reagrange and get the Shape Function: SF
SF=[s.c0,s.c1,s.c2,s.c3,s.c4,s.c5];
for kp=1:length(um)
    C=coeffs(un,um(kp));
    SF(kp)=C(2);   
end
SF2=diff(SF,x,2);
% calculate the stiffness matrix
Ke=int (EI*transpose(SF2)*SF2,x,0,L)
