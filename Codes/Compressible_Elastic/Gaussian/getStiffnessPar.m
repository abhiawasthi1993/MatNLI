function [Kmat,Mmat] = getStiffnessPar(mesh,Evec,Nu,Rho)

P   = mesh.P;
NCA = mesh.NCA;
nu = size(P,1);
ne = size(NCA,1);
ngp = 2;
[GP,W] = Gauss(ngp);

parfor e = 1:ne
    a = NCA(e,:);
    xe = P(a,1);
    ye = P(a,2);

    Ee = Evec(a,1);
    
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
        %x=N*xe;
        %y=N*ye;
              
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
        
        E = N*Ee;
        Coef = E/((1+Nu)*(1-2*Nu));
        CMat = Coef * [1-Nu Nu 0; 
                      Nu 1-Nu 0; 
                      0 0 (1-2*Nu)/2];
        RMat= Rho * [1 0;
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