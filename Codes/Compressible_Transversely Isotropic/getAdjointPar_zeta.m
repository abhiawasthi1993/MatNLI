function Gmat = getAdjointPar_zeta(X,mesh,A,U)

P = mesh.P;
NCA = mesh.NCA;
nu = size(P,1);
ne = size(NCA,1);
ngp = 2;
[GP,W] = Gauss(ngp);

Gvec = X(1:nu,1);
Zeta = X(nu+[1:nu],1);
Phi = X(2*nu+[1:nu],1);

Etheta = A(1,:);
T = getTransformMatrix2D(Etheta);

parfor e = 1:ne
    a = NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);
    locu = zeros(1,8);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;
    locmu = a;
    mu = Gvec(a,1);
    phi = Phi(locmu,1);
    zeta = Zeta(locmu,1);
    Ue = U(locu,1);

    Ge = zeros(8,4);

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
        
        G = N*mu;
        zetaE = N*zeta;
        phiE = N*phi;
        K = 0.5;

        c11 = (16/9)*G;
        c12 = -(8/9)*G;
        c22 = (4/9)*G;
        c44 = 0;

        Cmat = [c11 c12  0;
                c12 c22  0;
                 0   0  c44];

        CMat = T'*Cmat*T;
        
        ge = zeros(8,4);
        ge = B'*CMat*(B*Ue)*N;

        Ge = Ge + w*detJac*ge;
    end
    [i,j] = meshgrid(locu,locmu);
    I(e,:) = i(:);
    J(e,:) = j(:);
    G1(e,:) = reshape(Ge',[numel(Ge) 1]);
end

I = reshape(I, [numel(I) 1]);
J = reshape(J, [numel(J) 1]);
G1 = reshape(G1, [numel(G1) 1]);

Gmat = sparse(I,J,G1);