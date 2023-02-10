function [Mmat,Kmat,Gmat,Hmat] = getStiffnessPar(mesh,Gvec,Nu,Rho)

P   = mesh.P;
NCA = mesh.NCA;
ne = size(NCA,1);
ngp = 2;
[GP,W] = Gauss(ngp);

parfor e = 1:ne
    a = NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);
    Ge = Gvec(a,1);
    locu = zeros(1,8);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;
    locp = e;
    
	Me = zeros(8,8);
	Ke = zeros(8,8);
	Gme = zeros(8,1);
	He = zeros(1,1);

    for g = 1:size(GP,1)
        r = GP(g,1);
        s = GP(g,2);
        w = W(g);
        
        [N,DNr,DNs] = getQ4shp(r,s);
        x = N*xe;
        y = N*ye;
        
        E = 2*(N*Ge)*(1+Nu);
        
        Jac=[DNr; DNs]*[xe ye];
        detJac = det(Jac);
        DN = Jac\[DNr;DNs];
        
        Nd = zeros(2,8);
        Nd(1,1:2:8) = N(1,:);
        Nd(2,2:2:8) = N(1,:);
        
        B = zeros(3,8);
        B(1,1:2:8) = DN(1,:);
        B(2,2:2:8) = DN(2,:);
        B(3,1:2:8) = DN(2,:);
        B(3,2:2:8) = DN(1,:);
        
        mu = E/(2*(1+Nu));
        Kb = E/(3*(1-2*Nu));
        
        I0 = 0.5*[2 0 0;0 2 0; 0 0 1];
        m = [1 1 0]';
        Dd = 2*mu*(I0 - 1/3*m*m');

        me = Rho*Nd'*Nd;
        ke = B'*Dd*B;
        ge = B'*m*1;
        he = 1*1/Kb;
        
        Me = Me + w*detJac*me;
        Ke = Ke + w*detJac*ke;
        Gme = Gme + w*detJac*ge;
        He = He + w*detJac*he;
    end
	[i,j] = meshgrid(locu,locu);
	I(e,:) = i(:);
    J(e,:) = j(:);
    K(e,:) = locp;
    [x,y] = meshgrid(locu,locp);
    X(e,:) = x(:);
    Y(e,:) = y(:);
    
    M1(e,:) = reshape(Me,[numel(Me) 1]);
    K1(e,:) = reshape(Ke,[numel(Ke) 1]);
    G1(e,:) = reshape(Gme.',[numel(Gme) 1]);
    H1(e,:) = He;
end

I = reshape(I,[numel(I) 1]);
J = reshape(J,[numel(J) 1]);
X = reshape(X,[numel(X) 1]);
Y = reshape(Y,[numel(Y) 1]);

M1 = reshape(M1,[numel(M1) 1]);
K1 = reshape(K1,[numel(K1) 1]);
G1 = reshape(G1,[numel(G1) 1]);
Mmat = sparse(I,J,M1);
Kmat = sparse(I,J,K1);
Gmat = sparse(X,Y,G1);
Hmat = sparse(K,K,H1);
