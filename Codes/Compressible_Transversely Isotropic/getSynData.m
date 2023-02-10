%% COMPRESSIBLE TRANSVERSELY ISOTROPIC
clc
clear
close all

%% MATERIAL PARAMETERS
freq = 50;
Omega = 2*pi*freq;
Rho = 1000e-12;

%% GEOMETRY
lx = 120;
ly = 120;
nx = 51;
ny = 51;

%% MESHING
double = 0;                         % Put '1' for generating data with two BCs
mesh = getMesh(lx,ly,nx,ny,double);
nn = mesh.nn;
ne = mesh.ne;

dDoF = 2*nn;

%% MATERIAL PROPERTY VECTORS
Mu = getMu(mesh);
Zeta = getZeta(mesh);
Phi = getPhi(mesh);

Gvec = [Mu; Zeta; Phi];

plotSol(mesh.P,mesh.NCA,Mu,1,'\mu (MPa)');
plotSol(mesh.P,mesh.NCA,Phi,2,'\phi');
plotSol(mesh.P,mesh.NCA,Zeta,3,'\zeta');

%% FIBER ORIENTATION
beta = 30;                          % Angle (in degrees)
A = getFiberOrien(mesh, beta);      % Orthogonal basis corresponding to fiber orientation

%% SYSTEM MATRICES
[Kmat,Mmat] = getStiffnessPar(mesh,Gvec,A,Rho);
GK = Kmat - Omega^2*Mmat;
GF = zeros(dDoF,1);

%% BOUNDARY CONDITIONS
tdof = 1:dDoF;

switch(mesh.double)
    case 1
        kdof1 = sort([mesh.DBC1 mesh.hDBC1]);
        udof1 = setdiff(tdof,kdof1);

        kdof2 = sort([mesh.DBC2 mesh.hDBC2]);
        udof2 = setdiff(tdof,kdof2);

        Sol1 = zeros(2*nn,1);
        Sol1(mesh.DBC1,1) = 1;
        Sol1(mesh.hDBC1,1) = 0;

        Sol2 = zeros(2*nn,1);
        Sol2(mesh.DBC2,1) = 1;
        Sol2(mesh.hDBC2,1) = 0;

    otherwise
        kdof = sort([mesh.DBC mesh.hDBC]);
        udof = setdiff(tdof,kdof);

        Sol = zeros(2*nn,1);
        Sol(mesh.DBC,1) = 1;
        Sol(mesh.hDBC,1) = 0;
end

%% SOLVE FOR THE UNKNOWN DoF
switch(mesh.double)
    case 1
        Sol1(udof1,1) = GK(udof1,udof1)\(GF(udof1,1)-GK(udof1,kdof1)*Sol1(kdof1,1));
        Sol2(udof2,1) = GK(udof2,udof2)\(GF(udof2,1)-GK(udof2,kdof2)*Sol2(kdof2,1));
        Sol = [Sol1; Sol2];

    otherwise
        Sol(udof,1) = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*Sol(kdof,1));
end

%% POST-PROCESSING
switch(mesh.double)
    case 1
        Um = Sol;
        Ux1 = Um(1:2:dDoF);
        Uy1 = Um(2:2:dDoF);

        Ux2 = Um(dDoF+[1:2:dDoF]);
        Uy2 = Um(dDoF+[2:2:dDoF]);

    otherwise
        Um = Sol;
        Ux = Um(1:2:end);
        Uy = Um(2:2:end);
end

%% PLOT SOLUTION DATA
switch(mesh.double)
    case 1
        plotSol(mesh.P,mesh.NCA,Ux1,14,'X-Disp-TT');
        plotSol(mesh.P,mesh.NCA,Uy1,15,'Y-Disp-TT');
        plotSol(mesh.P,mesh.NCA,Ux2,16,'X-Disp-LTB');
        plotSol(mesh.P,mesh.NCA,Uy2,17,'Y-Disp-LTB');

    otherwise
        plotSol(mesh.P,mesh.NCA,Ux,14,'X-Disp-TT');
        plotSol(mesh.P,mesh.NCA,Uy,15,'Y-Disp-TT');
end

%% SAVING DATA
switch(double)
    case 1
        save('data_double.mat', 'mesh', 'Um', 'Mu', 'A', 'Omega', 'Rho', 'Gvec')
    otherwise
        save('data.mat', 'mesh', 'Um', 'Mu', 'A', 'Omega', 'Rho', 'Gvec')
end

%% SAVE SCALING DATA
D1 = (1/200)*eye(nn,nn);    % Corresponding to \mu
D2 = (1/0.5)*eye(nn,nn);    % Corresponding to \zeta
D3 = (1/5)*eye(nn,nn);      % Corresponding to \phi

ct1 = (1/100)*ones(nn,1);   % Corresponding to \mu
ct2 = 0.5*ones(nn,1);       % Corresponding to \zeta
ct3 = 0.2*ones(nn,1);       % Corresponding to \phi

D = sparse(blkdiag(D1, D2, D3));
Ct = [ct1; ct2; ct3];

save('ScData.mat','D','Ct')
