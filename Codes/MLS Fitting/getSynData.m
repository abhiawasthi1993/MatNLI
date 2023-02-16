%% COMPRESSIBLE ELASTIC
clc
close all
clear all

addpath('MLS')

%% MATERIAL PARAMETERS
Nu = 0.3;
Rho = 1000e-12;
freq = 80;
Omega = 2*pi*freq;

%% GEOMETRY
lx = 20;
ly = 20;
nx = 51;
ny = 51;

%% MESHING
mesh = getMesh(lx,ly,nx,ny);
nn = mesh.nn;
ne = mesh.ne;
dDoF = 2*nn;

%% MATERIAL PROPERTY VECTOR
Evec = getMatProp(mesh);
plotSol(mesh.P,mesh.NCA,Evec,1,'Young''s Modulus');

%% SYSTEM MATRICES
[Kmat,Mmat] = getStiffnessPar(mesh,Evec,Nu,Rho);

GK = Kmat - Omega^2*Mmat;
GF = zeros(dDoF,1);

%% BOUNDARY CONDITIONS
tDoF = 1:(2*nn);
kDoF = sort([mesh.DBC mesh.hDBC]);
uDoF = setdiff(tDoF,kDoF);

Sol = zeros(dDoF,1);
Sol(mesh.hDBC,1) = 0.0;
Sol(mesh.DBC,1) = 1;

%% SOLVE FOR THE UNKNOWN DoF
Sol(uDoF,1) = GK(uDoF,uDoF)\(GF(uDoF,1)-GK(uDoF,kDoF)*Sol(kDoF,1));

Ux = Sol(1:2:dDoF,1);
Uy = Sol(2:2:dDoF,1);

plotSol(mesh.P,mesh.NCA,Ux,12,'Clean X-Disp');
plotSol(mesh.P,mesh.NCA,Uy,13,'Clean Y-Disp');

%% ADDING WHITE GAUSSIAN NOISE
snr = 30;           % Signal to Noise ratio

Uxn = awgn(Ux,snr,'measured');    % Add noise in Ux
Uyn = awgn(Uy,snr,'measured');    % Add noise in Uy

Un = zeros(dDoF,1);  % Noisy displacement vector
Un(1:2:dDoF,1) = Uxn;
Un(1:2:dDoF,1) = Uyn;

plotSol(mesh.P,mesh.NCA,Uxn,22,'Noisy X-Disp');
plotSol(mesh.P,mesh.NCA,Uyn,23,'Noisy Y-Disp');

%% MLS FITTING
deltaX = 0.75;
deltaY = 0.75;

Um = zeros(dDoF,1);
Um(1:2:dDoF,1) = FitMLSSurface(mesh.P(:,1),mesh.P(:,2),Uxn,deltaX,deltaY);
Um(2:2:dDoF,1) = FitMLSSurface(mesh.P(:,1),mesh.P(:,2),Uyn,deltaX,deltaY);

plotSol(mesh.P,mesh.NCA,Um(1:2:dDoF,1),32,'Fitted X-Disp');
plotSol(mesh.P,mesh.NCA,Um(2:2:dDoF,1),33,'Fitted Y-Disp');

%% SAVING DATA
save('data.mat', 'mesh', 'Um', 'Evec', 'Omega', 'Nu', 'Rho')
