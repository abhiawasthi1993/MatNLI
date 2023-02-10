% Linearly elastic - Nearly Incompressible
clc
clear all
close all

%% MATERIAL PARAMETERS
freq = 40;
Omega = 2*pi*freq;
Nu = 0.4995;
Rho = 1000e-12;

%% GEOMETRY
lx = 120;
ly = 120;
nx = 101;
ny = 101;

%% MESHING
mesh = getMesh(lx,ly,nx,ny);
nn = mesh.nn;
ne = mesh.ne;

dDoF = 2*nn;
prDoF = ne;

%% MATERIAL PROPERTY VECTOR
Evec = getMatProp(mesh);
plotSol(mesh.P,mesh.NCA,Evec,1,'Young''s modulus (MPa)');

%% SYSTEM MATRICES
[M,K,G,H] = getStiffnessPar(mesh,Evec,Nu,Rho);
GK = [K-Omega^2*M G; G' -H];
GF = zeros(dDoF+prDoF,1);

%% BOUNDARY CONDITIONS
tDoF = 1:(dDoF + prDoF);
kDoF = unique([2*mesh.DBC-1 2*mesh.DBC]);
uDoF = setdiff(tDoF,kDoF);

Sol = zeros(dDoF+prDoF,1);
Sol(2*mesh.DBC,1) = 0;
Sol(2*mesh.DBC-1,1) = 1;

%% SOLVE FOR THE UNKNOWN DoF
Sol(uDoF,1) = GK(uDoF,uDoF)\(GF(uDoF,1) - GK(uDoF,kDoF)*Sol(kDoF,1));

%% POST-PROCESSING
Um = Sol(1:dDoF,1);
Pr = Sol([(dDoF+1):end],1);

Ux = Um(1:2:end);
Uy = Um(2:2:end);

%% PLOT SOLUTION DATA
plotSol(mesh.P,mesh.NCA,Ux,2,'X-Disp');
plotSol(mesh.P,mesh.NCA,Uy,3,'Y-Disp');
plotPress(mesh.P,mesh.NCA,Pr,4,'Pressure');

%% SAVING DATA
save('data.mat', 'mesh', 'Um', 'Pr', 'Evec', 'Nu', 'Rho', 'Omega')
