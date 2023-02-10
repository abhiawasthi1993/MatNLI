function Gmat = getAdjointPar(Gvec,mesh,Nu,U,Pr)

P = mesh.P;
NCA = mesh.NCA;
nu = size(P,1);
ne = size(NCA,1);
np = ne;
ngp = 2;
[GP,W] = Gauss(ngp);

parfor e = 1:ne
    a = NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);
    locu = zeros(1,8);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;
    locp = a;
    Ge = Gvec(locp,1);
    Ue = U(locu,1);
    Pre = Pr(e,1);

    Ge1 = zeros(8,4);
    Ge2 = zeros(1,4);

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
        
        E = 2*(N*Ge)*(1+Nu);
        mu = (N*Ge);
        coef = 1/(3*(1-2*Nu));
        Kb = E*coef;
		
		I0 = 0.5*[2 0 0;0 2 0; 0 0 1];
        m = [1 1 0]';
        Dd = 2*(I0 - 1/3*m*m');
		
		ge = zeros(8,4);
        ge = B'*Dd*(B*Ue)*N;
		he = zeros(1,4);
        he = 1*Pre*(1/Kb^2)*coef*N*2*(1+Nu);

        Ge1 = Ge1 + w*detJac*ge;
		Ge2 = Ge2 + w*detJac*he;
    end
    [i,j] = meshgrid(locu,locp);
    I(e,:) = i(:);
    J(e,:) = j(:);
    [x,y] = meshgrid(e,locp);
    X(e,:) = x(:);
    Y(e,:) = y(:);
    G1(e,:) = reshape(Ge1.',[numel(Ge1) 1]);
	G2(e,:) = reshape(Ge2.',[numel(Ge2) 1]);
end

I = reshape(I, [numel(I) 1]);
J = reshape(J, [numel(J) 1]);
X = reshape(X,[numel(X) 1]);
Y = reshape(Y,[numel(Y) 1]);

G1 = reshape(G1, [numel(G1) 1]);
G2 = reshape(G2, [numel(G2) 1]);

Gmat1 = sparse(I,J,G1);
Gmat2 = sparse(X,Y,G2);

Gmat = [Gmat1; Gmat2];