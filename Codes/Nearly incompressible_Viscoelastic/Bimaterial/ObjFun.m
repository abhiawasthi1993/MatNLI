function [Phi,Grad] = ObjFun(X,Um,mesh,Nu,Rho,Omega)

% DBC = mesh.BC;
DBC = unique([mesh.LB mesh.RB mesh.TB]);
nn  = mesh.nn;
np  = mesh.ne;

Gvec = X(1:nn,1) + 1i*X(nn+[1:nn],1);

%% BOUNDARY CONDITIONS
tdof = 1:(2*nn+np);			% Total degrees of freedom (dof)
kdof = [2*DBC-1 2*DBC];     % Known dof i.e bottom+left+right+top boundary dof
udof = setdiff(tdof,kdof);	% Unknown degrees of freedom

[M,K,C,V] = getStiffnessPar(mesh,Gvec,Nu,Rho);
GK = [K-Omega^2*M C; C' -V];
GF = zeros(2*nn+np,1);

Sol = zeros(2*nn+np,1);
Sol(kdof,1) = Um(kdof,1);
Sol(udof,1)  = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*Sol(kdof,1));

U = Sol(1:2*nn,1);
Pr = Sol((2*nn+1):end);

err = (U-Um);
Phi = 0.5*(err'*err);

if nargout > 1
    GF = zeros(2*nn+np,nn);
    GF(1:2*nn,1) = -conj(U-Um);
    GF((2*nn+1):end,1) = 0;
    
    W = zeros(2*nn+np,1);
    W(kdof,1) = 0.0;
    W(udof,1)  = GK(udof,udof)\(GF(udof,1));
    
    Gmat = getAdjointPar(Gvec,mesh,Nu,U,Pr);
    
    ADGrad(:,1) = conj(W.'*Gmat);
    Grad(1:nn,1) = real(ADGrad);
    Grad(nn+[1:nn],1) = imag(ADGrad);
end

