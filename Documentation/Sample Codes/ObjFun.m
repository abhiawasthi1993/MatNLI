function [Phi,AdGrad] = ObjFun(X,Um,mesh,Nu,Rho,Omega)

DBC = mesh.BC;          % Boundary nodes 
nn  = mesh.nn;          % Total number of nodes
np  = mesh.ne;          % Total number of elements

%% CONVERT REAL - VALUED STIFFNESS VECTOR TO COMPLEX - VALUED
Gvec = X(1:nn,1) + 1i*X(nn+[1:nn],1);

%% SYSTEM MATRICES
[M,K,C,V] = getStiffnessPar(mesh,Gvec,Nu,Rho);
GK = [K-Omega^2*M C; C' -V];
GF = zeros(2*nn+np,1);
Sol = zeros(2*nn+np,1);

%% BOUNDARY CONDITIONS
tdof = 1:(2*nn+np);	
kdof = unique([2*DBC-1 2*DBC]);
udof = setdiff(tdof,kdof);

%% SOLVE FOR UNKNOWNS
Sol(kdof,1) = Um(kdof,1);
Sol(udof,1) = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*Sol(kdof,1));

U = Sol(1:2*nn,1);
Pr = Sol((2*nn+1):end);

%% COMPUTE OBJECTIVE FUNCTION
err = (U-Um);
Phi = 0.5*(err'*err);

%% ADJOINT BASED GRADIENT COMPUTATION
if nargout > 1
    % FOrcing function for difference driven problem
    GF = zeros(2*nn+np,nn);
    GF(1:2*nn,1) = -conj(U-Um);
    GF((2*nn+1):end,1) = 0;
    
    % Computing unknown weights - Step 2
    W = zeros(2*nn+np,1);
    W(kdof,1) = 0.0;
    W(udof,1) = GK(udof,udof)\(GF(udof,1));
    
    % Assembly of [K'] matrix
    Gmat = getAdjointPar(Gvec,mesh,Nu,U,Pr);
    
    % Compute gradient  - Step 3
    Grad(:,1) = conj(W.'*Gmat);
    AdGrad(1:nn,1) = real(Grad);
    AdGrad(nn+[1:nn],1) = imag(Grad);
end

