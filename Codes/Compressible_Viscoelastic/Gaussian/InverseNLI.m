%% Compressible Viscoelastic
clc
close all
clear all

%% LOAD DATA
load data;
nn = mesh.nn;

%% INITIAL GUESS
X0(1:nn,1) = 0.01*ones(nn,1);          % Initial distribution of the Young's modulus as a vector
X0(nn+[1:nn],1) = 0.001*ones(nn,1);    % Initial distribution of the Young's modulus as a vector

%% OBJECTIVE FUNCTION DEFINITION
newfun = @(X) ObjFun(X,Um,mesh,Nu,Rho,Omega);

%% OPTIMIZATION
options = optimoptions('fminunc',...
    'PlotFcn',@optimplotfval,...
    'OutputFcn', @plotPar,...
    'Algorithm','quasi-newton',...
    'HessUpdate','bfgs',...
    'Display','iter-detailed',...
    'SpecifyObjectiveGradient',true,...
    'UseParallel',true,...
    'MaxIter',10000,...
    'MaxFunEvals',1e25,...
    'OptimalityTolerance',eps,...
    'TolX',eps,...
    'TolFun',eps);

%% CALL THE OPTIMIZER
tic
X = fminunc(newfun,X0,options);
toc

%%
plotSol(mesh.P,mesh.NCA,real(Gvec),3,'Predicted Storage Modulus');
plotSol(mesh.P,mesh.NCA,imag(Gvec),4,'Predicted Loss Modulus');

plotSol(mesh.P,mesh.NCA,X(1:mesh.nn,1),5,'Actual Storage Modulus');
plotSol(mesh.P,mesh.NCA,X(mesh.nn+[1:mesh.nn],1),6,'Actual Loss Modulus');
