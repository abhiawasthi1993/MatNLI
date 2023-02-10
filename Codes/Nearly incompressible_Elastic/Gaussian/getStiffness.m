function [M,K,G,H] = getStiffness(P,NCA,Evec,Nu,Rho)

nu = size(P,1);
ne = size(NCA,1);
np = ne;
ngp = 2;
[GP,W] = Gauss(ngp);

M = sparse(2*nu,2*nu);
K = sparse(2*nu,2*nu);
G = sparse(2*nu,np);
H = sparse(np,np);

for e=1:ne
    a=NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);
    Ee = Evec(a,1);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;
    locp = e;

    me_ser = zeros(8,8);
    ke_ser = zeros(8,8);
    ge_ser = zeros(8,1);
    he_ser = 0;

    for g=1:size(GP,1)
        r=GP(g,1);
        s=GP(g,2);
        w=W(g);
        
        [N,DNr,DNs] = getQ4shp(r,s);
        x=N*xe;
        y=N*ye;
        
        E = N*Ee;
        %E = GetYoungsMod(x,y)
        
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

        me_ser = me_ser + w*detJac*me;
        ke_ser = ke_ser + w*detJac*ke;
        ge_ser = ge_ser + w*detJac*ge;
        he_ser = he_ser + w*detJac*he;
        
        
        M(locu,locu)  = M(locu,locu)  + w*detJac*me;
        K(locu,locu)  = K(locu,locu)  + w*detJac*ke;
        G(locu,locp)  = G(locu,locp)  + w*detJac*ge;
        H(locp,locp)  = H(locp,locp)  + w*detJac*he;
    end
end
