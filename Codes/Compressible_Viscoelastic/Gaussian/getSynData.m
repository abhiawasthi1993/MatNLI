%% Compressible Viscoelastic
clc
clear all
close all

%% MATERIAL PARAMETERS
freq = 50;
Omega = 2*pi*freq;
Nu = 0.3;
Rho = 1000e-12;

%% GEOMETRY
lx = 120;
ly = 120;
nx = 51;
ny = 51;

%% MESHING & SOLUTION
mesh = getMesh(lx,ly,nx,ny);
nn = mesh.nn;
dDoF = 2*nn;

%% COMPLEX-VALUED MATERIAL PROPERTY VECTOR
Gvec = getMatProp(mesh);

plotSol(mesh.P,mesh.NCA,real(Gvec),1,'Storage modulus (MPa)');
plotSol(mesh.P,mesh.NCA,imag(Gvec),2,'Loss modulus (MPa)');

%% SYSTEM MATRICES
[Kmat,Mmat] = getStiffnessPar(mesh,Gvec,Nu,Rho);
GK = Kmat - Omega^2*Mmat;
GF = zeros(dDoF,1);

%% BOUNDARY CONDITIONS
tDoF = 1:dDoF;
kDoF = unique([2*mesh.DBC-1 2*mesh.DBC]);	   % Known degrees of freedom
uDoF = setdiff(tDoF,kDoF);

Sol = zeros(dDoF,1);
Sol(2*mesh.DBC,1) = 0;
Sol(2*mesh.DBC-1,1) = 1;

%% SOLVE FOR THE UNKNOWN DoF
Sol(uDoF,1)  = GK(uDoF,uDoF)\(GF(uDoF,1)-GK(uDoF,kDoF)*Sol(kDoF,1));

%% POST-PROCESSING
Um = Sol;

Ux = Um(1:2:end,1);
Uy = Um(2:2:end,1);

%% PLOTTING THE SOLUTION
plotSol(mesh.P,mesh.NCA,real(Ux),3,'Real X-Disp');
plotSol(mesh.P,mesh.NCA,real(Uy),4,'Real Y-Disp');
plotSol(mesh.P,mesh.NCA,imag(Ux),5,'Imag X-Disp');
plotSol(mesh.P,mesh.NCA,imag(Uy),6,'Imag Y-Disp');

%% SAVING DATA
save('data.mat', 'mesh', 'Gvec', 'Um', 'Omega', 'Nu', 'Rho')
