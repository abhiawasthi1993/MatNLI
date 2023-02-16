clc;
clear all;
close all;

%% MATERIAL PARAMETERS
Nu  = 0.4995;           % Poisson's ratio
Rho = 1000e-12;         % Density (kg/mm^-3)
freq = 100;             % Frequency (Hz)
Omega = 2*pi*freq;      % Angular frequency (rad)

%% GEOMETRY
lx = 120;               % Length of the domain along x-direction
ly = 120;               % Length of the domain along y-direction
nx = 51;                % No of nodes along x-direction
ny = 51;                % No of nodes along y-direction

%% MESHING - MESH AND BOUNDARY NODES INFORMATION
mesh = getMesh(lx,ly,nx,ny);
nn = mesh.nn;                   % Number of nodes
ne = mesh.ne;                   % Number of elements
dDoF = 2*nn;                    % Number of displacement DoF
prDoF = ne;                     % Number of pressure DoF

%% COMPLEX VALUED MATERIAL PROPERTY VECTOR
Gvec = GetMatProp(mesh);        % Complex shear modulus (MPa)

plotSol(mesh.P,mesh.NCA,real(Gvec),1,'Storage Shear Modulus')
plotSol(mesh.P,mesh.NCA,imag(Gvec),2,'Loss Shear Modulus')

%% COMPLEX VALUED SYSTEM MATRICES
[M,K,C,V] = getStiffnessPar(mesh,Gvec,Nu,Rho);
GK = [K-Omega^2*M C; C' -V];        % Global system matrix
GF = zeros(dDoF+prDoF,1);           % Global force vector

%% BOUNDARY CONDITIONS
tDoF = 1:(dDoF+prDoF);			        % Total Degrees of Freedom (dof)
kDoF = unique([mesh.DBC mesh.NBC]);		% Known DoF
uDoF = setdiff(tDoF,kDoF);	            % Unknown DoF

Sol = zeros(dDoF+prDoF,1);              % Initialize solution vector
Sol(mesh.DBC,1) = 0;                    % Homogeneous BC
Sol(mesh.NBC,1) = 1;                    % Dirichlet BC (mm)

%% SOLVE FOR THE UNKNOWN DoF
Sol(uDoF,1) = GK(uDoF,uDoF)\(GF(uDoF,1) - GK(uDoF,kDoF)*Sol(kDoF,1));

%% POST-PROCESSING
Um = Sol(1:dDoF,1);                 % Displacement solution
Pr = Sol([(dDoF+1):end],1);         % Pressure solution

Ux = Um(1:2:end);                   % Displacement in x-direction
Uy = Um(2:2:end);                   % Displacement in y-direction

%% Plot the data
plotSol(mesh.P,mesh.NCA,real(Ux),3,'Real - Ux');    % Plot x- displacement
plotSol(mesh.P,mesh.NCA,real(Uy),4,'Real - Uy');    % Plot y- displacement
plotPress(mesh.P,mesh.NCA,real(Pr),5,'Real - Pressure');          % Plot pressure

%% SAVING DATA
save('data.mat', 'mesh', 'Gvec', 'Um', 'Pr', 'Omega', 'Rho', 'Nu') 
