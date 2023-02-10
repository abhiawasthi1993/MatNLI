%% NLI - Compressible Elatic
clc
close all
clear

%% LOAD DATA
load data;
nn = mesh.nn;

%% INITIAL GUESS
X0 = 0.01*ones(nn,1);

%% OBJECTIVE FUNCTION DEFINITION
newfun = @(X) ObjFun(X,Um,mesh,Nu,Rho,Omega);

[~,GRAD] = ObjFun(X0,Um,mesh,Nu,Rho,Omega);

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
'MaxFunEvals',1e25,...
'OptimalityTolerance',eps,...
'TolX',eps,...
'TolFun',eps);

%% CALL THE OPTIMIZER
tic
X = fminunc(newfun,X0,options);
toc

%% POST PROCESSING
plotSol(mesh.P,mesh.NCA,Evec,4,'Actual Young''s Modulus'); caxis([0.01 0.02])
plotSol(mesh.P,mesh.NCA,X,5,'Predicted Young''s Modulus'); caxis([0.01 0.02])
