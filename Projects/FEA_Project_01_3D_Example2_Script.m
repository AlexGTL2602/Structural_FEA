%%
% Title : Truss Analysis in 3-Dimensions (3D)
% Author: Feng Zhou, Zhao Zi Jie, Nio Zi Feng, Alejandro Rodriguez, Victor
% Date : 15/11/2021
%==============================================================
% NOTATIONS                                           UNITS
%==============================================================
% E = Young's Modulus                               (N/m-sq)
% L = Length                                        (m)
% A = Area                                          (m-sq)
% force = force vector (F1x,F1y...etc)              (Pa)
% U = displacement vector (u1x,u1y,u1z..etc)        (m)
% element = Element from Node-to-Node (Element 1-2)
% nodeCoordinate = Coordinates Lists
% stiffness = Stiffness Matrix                      (N/m)
%==============================================================
%                          READ ME
%==============================================================
% There are 3 Inputs by User:
%
%      i. UserInput_1
%        - Input:
%          E,A, Coordinates, Nodes, Elements 
%          E.g: Element 1-2 (Meaning Bar Element from Node1 to Node2)
%               element = [1 2]
%
%               element = [1 2; 1 3; 1 4;]
%               (Meaning Bar Element from 1-2, 2-3 & 1-4)
%      
%     ii. UserInput_2
%        - Input:
%          Applied Forces
%          E.g: Applied 10kN downward force at Node1 Y-direction
%               force(2) = -10000
%
%    iii. UserInput_3
%        - Input:
%          PrescribedDOFs (Intialise Boundary Conditon)
% 
%==============================================================
%                      CODE STARTS BELOW
%==============================================================
%% Clear All Previous Data
clear all; 
clc; 
close all;

%% UserInput_1 - Input E,A, Coordinates, Nodes, Elements 
% Member Connectivity
element=[1 2; 1 3; 1 4; 1 5]; % Element 1 2 = element1-2 (Bar from Node 1 to 2)
nodeCoordinate=[ 4 4 3;  0 4 0;  0 4 6;   4 0 3;  8 -1 1]
numberElements=size(element,1); % Don't Edit this line
numberNodes=size(nodeCoordinate,1); % Don't Edit this line
GDof=3*numberNodes; % Global Degree of Freedom

% Member Properties
E_vec=210000*ones(1,length(element)); 
A_vec=100*ones(1,length(element));
EA = E_vec.*A_vec;


%% Initialize Displacement & Force Vector
D_vector = zeros(GDof,1); % displacement vector
F_vector = zeros(GDof,1);



%% UserInput_2 - Input Applied Forces / Loads 
% applied load at node 1 in y direction
F_vector(2) = -10000; % F1y


%% Generate Stiffness Matrix
% Generate Global Stiffness Matrix
[stiffness] = formStiffness3DTruss(GDof,numberElements,element,nodeCoordinate,E_vec,A_vec)


%% UserInput_3 - Define Zero-Value Displacement
prescribedDof=[4:GDof]';


%% Compute Displacements & Stresses
% Displacements
[D_vector,F_vector]=solution(prescribedDof,stiffness,D_vector,F_vector)


outputDisplacementsReactions(D_vector,stiffness,GDof,prescribedDof);

% Stresses at elements
stresses3Dtruss(numberElements,element,nodeCoordinate,D_vector,E_vec)