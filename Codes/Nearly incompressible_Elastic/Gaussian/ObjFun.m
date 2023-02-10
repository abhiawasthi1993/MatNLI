function [Phi,ADGrad] = ObjFun(X,Um,mesh,Nu,Rho,Omega)

DBC = mesh.BC;
nn  = mesh.nn;
np  = mesh.ne;

%% SYSTEM MATRICES
[M,K,C,V] = getStiffnessPar(mesh,X,Nu,Rho);
GK = [K-Omega^2*M C; C' -V];
GF = zeros(2*nn+np,1);
Sol = zeros(2*nn+np,1);

%% BOUNDARY CONDITIONS
tDoF = 1:(2*nn+np);
kDoF = unique([2*DBC-1 2*DBC]);
uDoF = setdiff(tDoF,kDoF);

%% SOLVE FOR UNKNOWNS
Sol(kDoF,1) = Um(kDoF,1);
Sol(uDoF,1) = GK(uDoF,uDoF)\(GF(uDoF,1)-GK(uDoF,kDoF)*Sol(kDoF,1));

U = Sol(1:2*nn,1);
Pr = Sol((2*nn+1):end);

%% OBJECTIVE FUNCTION EVALUATION
err = (U - Um);
Phi = 0.5*(err'*err);

%% ADJOINT BASED GRADIENT COMPUTATION
if nargout > 1
    GF = zeros(2*nn+np,nn);
    GF(1:2*nn,1) = -(U-Um);
    GF((2*nn+1):end,1) = 0;
    
    W = zeros(2*nn+np,1);
    W(kDoF,1) = 0.0;
    W(uDoF,1)  = GK(uDoF,uDoF)\(GF(uDoF,1));
    
    Gmat = getAdjointPar(X,mesh,Nu,U,Pr);
    
    ADGrad(:,1) = W'*Gmat;
end

