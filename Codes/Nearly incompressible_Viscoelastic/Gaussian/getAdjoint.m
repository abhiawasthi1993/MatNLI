function Gmat = getAdjoint(Gvec,mesh,Nu,Rho,U,Pr)

P   = mesh.P;
NCA = mesh.NCA;
nu = size(P,1);
np = size(NCA,1);
ne = size(NCA,1);
ngp = 2;
[GP,W] = Gauss(ngp);

Gmat1 = sparse(2*nu,nu);
Gmat2 = sparse(ne,nu);

for e=1:ne
    a=NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);
    locu = zeros(1,8);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;
    locp = a;
    Ge = Gvec(locp,1);
    Ue = U(locu,1);
    Pre = Pr(e,1);
    for g=1:size(GP,1)
        r=GP(g,1);
        s=GP(g,2);
        w=W(g);
        
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
        he = 1*Pre*(1/Kb^2)*coef*N;

        Gmat1(locu,locp)  = Gmat1(locu,locp)   + w*detJac*ge; 
        Gmat2(e,locp)     = Gmat2(e,locp)      + w*detJac*he; 
    end
end
Gmat = [Gmat1; Gmat2];
