% NLI - Viscoelastic Nearly-Incompressible
clc
close all
clear all

%% LOAD DATA
load data
nn = mesh.nn;

%% INITIAL GUESS
X0(1:nn,1) = 0.01*ones(nn,1);       
X0(nn+[1:nn],1) = 0.001*ones(nn,1);

%% OBJECTIVE FUNCTION DEFINITION
newfun = @(X) ObjFun(X,Um,mesh,Nu,Rho,Omega);

%% SETTING OPTIONS FOR OPTIMIZAION MODULE
options = optimoptions('fminunc',...
'PlotFcn',@optimplotfval,...
'OutputFcn', @plotPar,...
'Algorithm','quasi-newton',...
'HessUpdate','bfgs',...
'Display','iter-detailed',...
'SpecifyObjectiveGradient',true,...
'UseParallel',true,...
'MaxIter',10000,...
'MaxFunEvals',1e50,...
'OptimalityTolerance',eps,...
'TolX',eps,...
'TolFun',eps);

%% CALL THE OPTIMIZER
tic
X = fminunc(newfun,X0,options);
toc

%% POST PROCESSING
plotsol(mesh.P,mesh.NCA,X(1:mesh.nn,1),5,'Predicted Storage Modulus');
plotsol(mesh.P,mesh.NCA,X(mesh.nn+[1:mesh.nn],1),6,'Predicted Loss Modulus');

plotsol(mesh.P,mesh.NCA,real(Gvec),7,'Actual Storage Modulus');
plotsol(mesh.P,mesh.NCA,imag(Gvec),8,'Actual Loss Modulus');

