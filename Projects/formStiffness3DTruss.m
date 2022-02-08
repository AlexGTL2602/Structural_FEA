%%
% Title : Assemblage / Assembly Function Code (3D)
% Author: Feng Zhou, Zhao Zi Jie, Nio Zi Feng, Alejandro Rodriguez, Victor
% Date : 15/11/2021
%==============================================================
% NOTATIONS                                           UNITS
%==============================================================
% E = Young's Modulus                               (N/m-sq)
% L = Length                                        (m)
% A = Area                                          (m-sq)
% theta = angle                                     (radian)
% force = force vector (F1x,F1y...etc)              (Pa)
% U = displacement vector (u1x,u1y,u1z..etc)        (m)
% length_element = Bar Length / Element Length      (m
% Cx = Cos(theta_x) 
% Cy = Cos(theta_y)
% Cz = Cos(theta_z)
% stiffness = Stiffness Matrix                      (N/m)
%==============================================================
%                      CODE STARTS BELOW
%==============================================================

%%
function [stiffness]=formStiffness3DTruss(GDof,numberElements,elementNodes,nodeCoordinate,E,A)

% Initialise the Global Stiffness Matrix with all Zeros
stiffness=zeros(GDof);

for i=1 : numberElements
    % elementDof: element degrees of freedom (Dof)
    indice=elementNodes(i,:); 


    elementDof=...
        [ indice(1)*3-2 indice(1)*3-1 indice(1)*3 indice(2)*3-2 indice(2)*3-1 indice(2)*3];
    
    % Note: Indice 1 = First Node X,Y,Z , Indice(2) = Second Node X,Y,Z 
    % Eg. If it's Element 1 - 2 ( Indice(1) will gives u x1,y1,z1,
    % Indice(2) will gives u x2,y2,z2)
    
    
    % Indice(1) gives us X, Y, Z
    % Indice(1),1 = X1
    % Indice(1),2 = Y1
    % Indice(1),3 = Z1

    % Indice(2) gives us X, Y, Z
    % Indice(2),1 = X2
    % Indice(2),2 = Y2
    % Indice(2),3 = Z2

    x1=nodeCoordinate(indice(1),1);
    y1=nodeCoordinate(indice(1),2);
    z1=nodeCoordinate(indice(1),3);

    x2=nodeCoordinate(indice(2),1);
    y2=nodeCoordinate(indice(2),2);
    z2=nodeCoordinate(indice(2),3);


    length_element=sqrt((x2-x1)^2 +  (y2-y1)^2  +  (z2-z1)^2 );

    % Cosine theta for x,y,z (This is the formulae, no angle, just taking
    % the x2-x1/Length or y2-y1/Length or z2-z1/Length
    Cx = (x2-x1)/length_element; % Lambda x

    Cy = (y2-y1)/length_element; % Lambda y

    Cz = (z2-z1)/length_element; % Lambda z

    % Transformation Matrix (The C*C CS stuff)
    T=[Cx*Cx Cx*Cy Cx*Cz; Cx*Cy Cy*Cy Cy*Cz; Cx*Cz Cy*Cz Cz*Cz]; 


    % Forming Stiffness matrix k1 (Local Stiffness Matrix)
    k1 = E(i)*A(i)/length_element.*[T -T; -T T]   %EA/L * [1 -1; 1 1] * T
    
   
    
    % Assemblage (Store the Local Stiffness Matrix into Global Stiffness Matrix)
    stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+k1;
    

end
