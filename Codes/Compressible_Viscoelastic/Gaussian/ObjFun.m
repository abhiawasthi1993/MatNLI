function [Phi,Grad] = ObjFun(X,Um,mesh,Nu,Rho,Omega)

DBC = mesh.BC;
nn  = mesh.nn;
Gvec = X(1:nn,1) + 1i*X(nn+[1:nn],1);

%% BOUNDARY CONDITIONS
tdof = 1:(2*nn);        			% Total degrees of freedom (dof)
kdof = unique([2*DBC-1 2*DBC]);		% Known dof i.e bottom+left+top boundary dof
udof = setdiff(tdof,kdof);          % Unknown degrees of freedom

[Kmat, Mmat] = getStiffnessPar(mesh,Gvec,Nu,Rho);
GK = Kmat - Omega^2*Mmat;
GF = zeros(2*nn,1);

Sol = zeros(2*nn,1);
Sol(kdof,1) = Um(kdof,1);
Sol(udof,1)  = GK(udof,udof)\(GF(udof,1) - GK(udof,kdof)*Sol(kdof,1));

U = Sol(1:2*nn,1);

err = (U(udof)-Um(udof));
Phi = 0.5*(err'*err);

if nargout > 1
    GF = -conj(U-Um);
    
    W = zeros(2*nn,1);
    W(kdof,1) = 0.0;
    W(udof,1) = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*W(kdof,1));
    
    Gmat = getAdjointPar(Gvec,mesh,Nu,U);
    
    ADGrad(:,1) = conj(W.'*Gmat);
    Grad(1:nn,1) = real(ADGrad);
    Grad(nn+[1:nn],1) = imag(ADGrad);

end

