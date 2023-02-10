function Gmat = getAdjoint(Gvec,mesh,Nu,Rho,U)

P   = mesh.P;
NCA = mesh.NCA;
nu = size(P,1);
ne = size(NCA,1);
ngp = 2;
[GP,W] = Gauss(ngp);

Gmat = sparse(2*nu,nu);

for e = 1:ne
    a = NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);
    locu = zeros(1,8);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;
    locp = a;
    Ue = U(locu,1);

    for g = 1:size(GP,1)
        r = GP(g,1);
        s = GP(g,2);
        w = W(g);
        
        [N,DNr,DNs] = getQ4shp(r,s);
              
        Jac=[DNr; DNs]*[xe ye];
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
        
        Coef = 2*(1+Nu)/((1+Nu)*(1-2*Nu));
        CMat= Coef * [1-Nu Nu 0; 
                      Nu 1-Nu 0; 
                      0 0 (1-2*Nu)/2];
        
        ge = zeros(8,4);
        ge = B'*CMat*(B*Ue)*N;

        Gmat(locu,locp)  = Gmat(locu,locp) + w*detJac*ge;
    end
end
