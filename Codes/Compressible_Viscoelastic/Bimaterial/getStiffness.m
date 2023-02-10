function [Kmat, Mmat] = getStiffness(mesh,Gvec,Nu,Rho)

P = mesh.P;
NCA = mesh.NCA;

nu = size(P,1);
ne = size(NCA,1);
ngp = 2;
[GP,W] = Gauss(ngp);

Mmat = sparse(2*nu,2*nu);
Kmat = sparse(2*nu,2*nu);

for e = 1:ne
    a = NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);
    Ge = Gvec(a,1);
    locu = zeros(1,8);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;

    for g = 1 : size(GP,1)
        r = GP(g,1);
        s = GP(g,2);
        w = W(g);

        [N,DNr,DNs] = getQ4shp(r,s);

        E = 2*(N*Ge)*(1+Nu);

        Jac = [DNr; DNs]*[xe ye];
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

        Coef = E/((1+Nu)*(1-2*Nu));
        CMat = Coef * [1-Nu Nu 0;
                       Nu 1-Nu 0;
                       0 0 (1-2*Nu)/2];
        RMat = Rho * [1 0;
                      0 1];

        me = zeros(8,8);
        ke = zeros(8,8);
        me = Nd'*RMat*Nd;
        ke = B'*CMat*B;

        Mmat(locu,locu)  = Mmat(locu,locu)  + w*detJac*me;
        Kmat(locu,locu)  = Kmat(locu,locu)  + w*detJac*ke;
    end
end
