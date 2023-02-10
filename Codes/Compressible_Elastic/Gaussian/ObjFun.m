function [Phi,GradPhi] = ObjFun(X,Um,mesh,Nu,Rho,Omega)

DBC = mesh.BC;
nn  = mesh.nn;

%% BOUNDARY CONDITIONS
tdof = 1:(2*nn);
kdof = sort([2*DBC-1 2*DBC]);
udof = setdiff(tdof,kdof);

[Kmat,Mmat] = getStiffnessPar(mesh,X,Nu,Rho);
GK = Kmat-Omega^2*Mmat;
GF = zeros(2*nn,1);

U(kdof,1) = Um(kdof,1);
U(udof,1)  = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*U(kdof,1));

err = (U-Um);
Phi = 0.5*(err'*err);

if nargout > 1
    GF = -(U-Um);
    Gmat = getAdjointPar(X,mesh,Nu,U);
    W(kdof,1) = 0.0;
    W(udof,1) = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*W(kdof,1));
    GradPhi(:,1) = W'*Gmat;
end

