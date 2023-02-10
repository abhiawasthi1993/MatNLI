%% Compressible Transversely Isotropic
clc
close all
clear all

%% LOAD DATA
load data.mat;
nn = mesh.nn;

%% INITIAL GUESS FOR PARAMETERS (MU, ZETA, PHI)
X0 = zeros(3*nn,1);

%% OBJECTIVE FUNCTION DEFINITION
newfun = @(X) ObjFun(X,Um,mesh,A,Rho,Omega);

%% SETTING OPTIONS FOR OPTIMIZAION MODULE
options = optimoptions('fminunc',...
                        'PlotFcn',@optimplotfval,...
                        'OutputFcn', @plotPar,...
                        'Algorithm','quasi-newton',...
                        'HessUpdate','bfgs',...
                        'Display','iter-detailed',...
                        'SpecifyObjectiveGradient',true,...
                        'UseParallel',true,...
                        'MaxIter',20000,...
                        'MaxFunEvals',1e25,...
                        'OptimalityTolerance',1e-6,...
                        'TolX',eps,...
                        'TolFun',eps);

%% CALL THE OPTIMIZER
tic
X = fminunc(newfun,X0,options);
toc

%% POST-PROCESSING
load('ScData.mat','D','Ct')
X1 = D*X + Ct;

plotSol(mesh.P,mesh.NCA,X1(1:nn,1),3,'Predicted - \mu'); caxis([0.01; 0.015])
plotSol(mesh.P,mesh.NCA,X1(nn+[1:nn],1),4,'Predicted - \zeta'); caxis([0.5; 2.5])
plotSol(mesh.P,mesh.NCA,X1(2*nn+[1:nn],1),5,'Predicted - \phi'); %caxis([0.25; 0.4])

plotSol(mesh.P,mesh.NCA,Gvec(1:nn,1),6,'Actual - \mu'); caxis([0.01; 0.015])
plotSol(mesh.P,mesh.NCA,Gvec(nn+[1:nn],1),7,'Actual - \zeta'); caxis([0.5; 2.5])
plotSol(mesh.P,mesh.NCA,Gvec(2*nn+[1:nn],1),8,'Actual - \phi'); caxis([0.2; 0.4])
