% Nearly Incompressible
clc;
clear all;
close all;

%% MATERIAL PARAMETERS
freq    = 100;
Omega   = 2*pi*freq;
Nu      = 0.49;
Rho     = 1000e-12;

%% GEOMETRY
lx = 120;
ly = 120;
nx = 51;
ny = 51;

%% MESHING
mesh = getMesh(lx,ly,nx,ny);
nn = mesh.nn;
ne = mesh.ne;
dDoF = 2*nn;
prDoF = ne;

%% COMPLEX-VALUED MATERIAL PROPERTY VECTOR
Gvec = GetMatMod(mesh);

plotSol(mesh.P,mesh.NCA,real(Gvec),1,'Storage Shear Modulus (MPa)')
plotSol(mesh.P,mesh.NCA,imag(Gvec),2,'Loss Shear Modulus (MPa)')

%% COMPLEX VALUED SYSTEM MATRICES
[M,K,C,V] = getStiffnessPar(mesh,Gvec,Nu,Rho);
GK = [K-Omega^2*M C; C' -V];
GF = zeros(dDoF+prDoF,1);

%% BOUNDARY CONDITIONS
tDoF = 1:(dDoF+prDoF);
kDoF = [2*mesh.DBC-1 2*mesh.DBC];
uDoF = setdiff(tDoF,kDoF);

Sol = zeros(dDoF+prDoF,1);
Sol(2*mesh.DBC,1) = 0;
Sol(2*mesh.DBC-1,1) = 1;

%% SOLVE FOR THE UNKNOWNS
Sol(uDoF,1) = GK(uDoF,uDoF)\(GF(uDoF,1)-GK(uDoF,kDoF)*Sol(kDoF,1));

%% POST-PROCESSING
Um = Sol(1:dDoF,1);
Pr = Sol([(dDoF+1):end],1);

Ux = Um(1:2:end,1);
Uy = Um(2:2:end,1);

%% PLOT THE DATA
plotSol(mesh.P,mesh.NCA,real(Ux),3,'Real X-Disp');
plotSol(mesh.P,mesh.NCA,imag(Ux),4,'Imag X-Disp');
plotSol(mesh.P,mesh.NCA,real(Uy),5,'Real Y-Disp');
plotSol(mesh.P,mesh.NCA,imag(Uy),6,'Imag Y-Disp');
plotPress(mesh.P,mesh.NCA,real(Pr),7,'Real Pressure');
plotPress(mesh.P,mesh.NCA,imag(Pr),8,'Imag Pressure');

%% SAVING DATA
save('data.mat', 'mesh', 'Um', 'Pr', 'Gvec', 'Nu', 'Rho', 'Omega')
