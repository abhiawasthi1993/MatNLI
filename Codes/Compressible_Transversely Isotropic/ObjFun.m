function [ObjFun,GradPhi] = ObjFun(X1,Um,mesh,A,Rho,Omega)

BC = (mesh.BC)';
nn  = mesh.nn;

load('ScData.mat','D','Ct')

X = D*X1 + Ct;

[Kmat,Mmat] = getStiffnessPar(mesh,X,A,Rho);
GK = Kmat - Omega^2*Mmat;
GF = zeros(2*nn,1);

switch(mesh.double)
    case 1
        tdof = 1:(2*nn);
        kdof = sort([2*BC-1; 2*BC]);
        udof = setdiff(tdof,kdof);

        U1 = zeros(2*nn,1);
        U1(kdof,1) = Um(kdof,1);
        U1(udof,1) = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*U1(kdof,1));

        U2 = zeros(2*nn,1);
        U2(kdof,1) = Um(2*nn+[kdof],1);
        U2(udof,1) = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*U2(kdof,1));

        U = [U1; U2];

        err = (U - Um);
        ObjFun = 0.5*(err'*err);

        if nargout > 1
            Gmat1_tb = getAdjointPar_mu(X,mesh,A,U(1:2*nn,1));
            Gmat2_tb = getAdjointPar_zeta(X,mesh,A,U(1:2*nn,1));
            Gmat3_tb = getAdjointPar_phi(X,mesh,A,U(1:2*nn,1));

            Gmat1_lr = getAdjointPar_mu(X,mesh,A,U(2*nn+[1:2*nn],1));
            Gmat2_lr = getAdjointPar_zeta(X,mesh,A,U(2*nn+[1:2*nn],1));
            Gmat3_lr = getAdjointPar_phi(X,mesh,A,U(2*nn+[1:2*nn],1));

            Gmat_tb = [Gmat1_tb Gmat2_tb Gmat3_tb];
            Gmat_lr = [Gmat1_lr Gmat2_lr Gmat3_lr];

            Gmat = [Gmat_tb; Gmat_lr];

            GF = -(U-Um);

            W1 = zeros(2*nn,1);
            W1(kdof,1) = 0.0;
            GF1 = GF(1:2*nn,1);
            W1(udof,1) = GK(udof,udof)\(GF1(udof,1) - GK(udof,kdof)*W1(kdof,1));

            W2 = zeros(2*nn,1);
            W2(kdof,1) = 0.0;
            GF2 = GF(2*nn+[1:2*nn],1);
            W2(udof,1) = GK(udof,udof)\(GF2(udof,1) - GK(udof,kdof)*W2(kdof,1));

            W = [W1; W2];

            GradPhi1(:,1) = W'*Gmat;

            grad_mu = GradPhi1(1:nn,1);
            grad_chi = GradPhi1(nn+[1:nn],1);
            grad_phi = GradPhi1(2*nn+[1:nn],1);

            GradPhi1 = [grad_mu; grad_chi; grad_phi];
            GradPhi = D*GradPhi1;
        end
    otherwise
        tdof = 1:(2*nn);
        kdof = sort([2*BC-1 2*BC]);
        udof = setdiff(tdof,kdof);

        U(kdof,1) = Um(kdof,1);
        U(udof,1) = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*U(kdof,1));

        err = (U - Um);
        ObjFun = 0.5*(err'*err);

        if nargout > 1
            GF = -(U-Um);
            Gmat1 = getAdjointPar_mu(X,mesh,A,U);
            Gmat2 = getAdjointPar_zeta(X,mesh,A,U);
            Gmat3 = getAdjointPar_phi(X,mesh,A,U);

            Gmat = [Gmat1 Gmat2 Gmat3];

            W(kdof,1) = 0.0;
            W(udof,1) = GK(udof,udof)\(GF(udof,1)-GK(udof,kdof)*W(kdof,1));
            GradPhi1(:,1) = W'*Gmat;

            grad_mu = GradPhi1(1:nn,1);
            grad_chi = GradPhi1(nn+[1:nn],1);
            grad_phi = GradPhi1(2*nn+[1:nn],1);

            GradPhi1 = [grad_mu; grad_chi; grad_phi];
            GradPhi = D*GradPhi1;
        end
end

