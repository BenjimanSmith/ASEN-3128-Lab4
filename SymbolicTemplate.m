%% put given values here
m = .068; % mass (kg)
r = 0.06; % radius (m)
Ix = 6.8e-5; % Moment of inertia abt x (kgm^2)
Iy = 9.2e-5; % Moment of inertia abt y (kgm^2)
Iz = 1.35e-4; % Moment of inertia abt z (kgm^2)
k = .0024; % Nm/N
g=9.81; % gravity (m/s^2)
alpha = 2e-6;
nu = 1e-3;
syms F1 F2 F3 F4 % Initialize unknown variables
%% four equations in this case, four unknowms
eqn1 = (r/sqrt(2))*(F1 + F2 - F3 -F4) == 0; 
eqn2 = (r/sqrt(2))*(-F1 + F2 + F3 -F4) == 0;
eqn3 = k*(F1 - F2 + F3 -F4) == 0;
eqn4 = F1 + F2 + F3 + F4 == g*m;

[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], [F1, F2, F3, F4]); % puts both the equations and the unknown values in a matrix

SolutionVect = linsolve(A,B); % vector for solutions
SolutionVect = double(SolutionVect); % convert the symbolic solution vector to a vector of doubles
F1 = SolutionVect(1); % 1st unknown
F2 = SolutionVect(2); % 2nd unknown
F3 = SolutionVect(3); % 3rd unknown
F4 = SolutionVect(4); % 4th unknown