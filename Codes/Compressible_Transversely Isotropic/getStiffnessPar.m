function [Kmat,Mmat] = getStiffnessPar(mesh,X,A,Rho)

P   = mesh.P;
NCA = mesh.NCA;
nu = size(P,1);
ne = size(NCA,1);
ngp = 2;
[GP,W] = Gauss(ngp);

Gvec = X(1:nu,1);
Chi = X(nu+[1:nu],1);
Phi = X(2*nu+[1:nu],1);

% Constatnt fiber orientation at every element (otherwise use lines 68-69)
Etheta = A(1,:);
T = getTransformMatrix2D(Etheta);

parfor e = 1:ne
    a = NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);
    Ge = Gvec(a,1);
    phi = Phi(a,1);
    chi = Chi(a,1);

    locu = zeros(1,8);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;

    Me = zeros(8,8);
    Ke = zeros(8,8);

    for g = 1:size(GP,1)
        r = GP(g,1);
        s = GP(g,2);
        w = W(g);

        [N,DNr,DNs] = getQ4shp(r,s);

        Jac = [DNr; DNs]*[xe ye];
        detJac = det(Jac);
        DN = Jac\[DNr;DNs];

        Nd = zeros(2,8);
        Nd(1,1:2:8) = N;
        Nd(2,2:2:8) = N;

        B = zeros(3,8);
        B(1,1:2:8) = DN(1,:);
        B(2,2:2:8) = DN(2,:);
        B(3,1:2:8) = DN(2,:);
        B(3,2:2:8) = DN(1,:);

        G = N*Ge;
        chiE = N*chi;
        phiE = N*phi;
        K = 0.5;            % in MPa

        c11 = K + (4/3)*G*(1+(4/3)*chiE);
        c12 = K - (2/3)*G*(1+(4/3)*chiE);
        c22 = K + (4/3)*G*(1+(1/3)*chiE);
        c44 = G*(1+phiE);

        Cmat = [c11 c12 0;
                c12 c22 0;
                 0   0  c44];

        % Etheta = A(e,:);
        % T = getTransformMatrix2D(Etheta);

        CMat = T'*Cmat*T;
        
        RMat = Rho*[1 0;
                    0 1];

        me = zeros(8,8);
        ke = zeros(8,8);
        me = Nd'*RMat*Nd;
        ke = B'*CMat*B;

        Me  = Me + w*detJac*me;
        Ke  = Ke + w*detJac*ke;
    end
    [i,j] = meshgrid(locu,locu);
    I(e,:) = i(:);
    J(e,:) = j(:);
    M1(e,:) = reshape(Me,[numel(Me) 1])';
    K1(e,:) = reshape(Ke,[numel(Ke) 1])';
end

I = reshape(I,[numel(I) 1]);
J = reshape(J,[numel(J) 1]);
M1 = reshape(M1,[numel(M1) 1]);
K1 = reshape(K1,[numel(K1) 1]);
Mmat = sparse(I,J,M1);
Kmat = sparse(I,J,K1);