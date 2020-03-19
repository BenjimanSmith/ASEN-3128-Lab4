% This script finds the values for derrivative gain controlling both
% the longitudional and lateral rates of the drone
% 
%
%   Author: Benjiman Smith
%   Collaborators: E. Owen, I. Quezada
%   Date: 1/25/2020
%
%% lateral control
clear all
close all
Lambda1 = -2; % 1st lambda value, found from lab calculations
Lambda2 = -12;% 2nd lambda value is 10x Lambda 1, in order to make lambda 1 dominant
Ix = 6.8e-5;   % moment of inertia in the x direction [kg m^2]
Iy = 9.2e-5; % Moment of inertia abt y (kgm^2)

syms K1 K2 % symbollically solving for K1 and K2
eqn1 = Lambda1^2 +Lambda1*(K1/Ix) + (K2/Ix) == 0; % 1st eigenvalue solution
eqn2 = Lambda2^2 +Lambda2*(K1/Ix) + (K2/Ix) == 0; % 2nd eigenvalue solution


[A,B] = equationsToMatrix([eqn1, eqn2], [K1, K2]); % convert both the equations and the unknowns into independent vectors

Solution = linsolve(A,B); % solve the system of equations
Solution = double(Solution); % convert symbolic solution vector to vector of doubles
K1 = Solution(1)/Ix; % bank control
K2 = Solution(2)/Ix; % bank control


%% Longitudional Control
syms K3 K4 % symbollically solving for K3 and K4
eqn3 = Lambda1^2 +Lambda1*(K3/Iy) + (K4/Iy) == 0; % 1st eigenvalue solution
eqn4 = Lambda2^2 +Lambda2*(K3/Iy) + (K4/Iy) == 0; % 1st eigenvalue solution


[C,D] = equationsToMatrix([eqn3, eqn4], [K3, K4]);

ForceVect2 = linsolve(C,D);
ForceVect2 = double(ForceVect2);
K3 = ForceVect2(1)/Iy; % elevation control
K4 = ForceVect2(2)/Iy; % elevation control