function [Mmat,Kmat,Cmat,Vmat] = getStiffnessPar(mesh,Gvec,Nu,Rho)

P   = mesh.P;                 % Nodal coordinates
NCA = mesh.NCA;               % Nodal Connectivity Array (NCA)
nu  = size(P,1);              % Number of nodes
ne  = size(NCA,1);            % Number of elements
np  = ne;                     % Number of Pressure DoFs

ngp = 2;                      % Number of Gauss points
[GP, W] = Gauss(ngp) ;  % Gauss point locations and weights

% Element level loop starts here
parfor e = 1:ne
    a = NCA(e,:);               % Extracting nodes for element 'e'
    xe = P(a,1);                % x-Coordinates of nodes of element 'e'
    ye = P(a,2);                % y-Coordinates of nodes of element 'e'
    Ge = Gvec(a,1);             % Nodal values of shear modulus for element 'e'
    locu = zeros(1,8);          % Location vector for global assembly
    locu(1:2:8) = 2*a-1;        % Locations for ux
    locu(2:2:8) = 2*a;          % Locations for uy
    locp = e;                   % Locations for pressure (p)
    
	Me = zeros(8,8);            % Element level mass matrix
	Ke = zeros(8,8);            % Element level stiffness matrix
	Ce = zeros(8,1);            % Mixed terms based on displacment and pressure
	Ve = zeros(1,1);            % Pressure dependent term

    % Computations at Gauss point level
    for g = 1:size(GP,1)
        r = GP(g,1);
        s = GP(g,2);
        w = W(g);
        
        % Get shape function and their derivatives
        [N,DNr,DNs] = getQ4shp(r,s);
        Jac = [DNr; DNs]*[xe ye];
        detJac = det(Jac);
        DN = Jac\[DNr;DNs];
        
        Nd = zeros(2,8);
        Nd(1,1:2:8) = N(1,:);
        Nd(2,2:2:8) = N(1,:);
        
        % Compute strain - displacement matrix (B)
        B = zeros(3,8);
        B(1,1:2:8) = DN(1,:);
        B(2,2:2:8) = DN(2,:);
        B(3,1:2:8) = DN(2,:);
        B(3,2:2:8) = DN(1,:);
        
        % Compute strain - displacement matrix (B)
        E = 2*(N*Ge)*(1+Nu);            % Young's modulus
        mu = E/(2*(1+Nu));              % Shear modulus
        Kb = E/(3*(1-2*Nu));            % Bulk modulus
        
        % Deviatoric part of constitutive matrix (Dd)
        I0 = 0.5*[2 0 0;0 2 0; 0 0 1];
        m = [1 1 0]';
        Dd = 2*mu*(I0 - 1/3*m*m');

        me = Rho*Nd'*Nd;
        ke = B'*Dd*B;
        ce = B'*m*1;
        ve = 1*1/Kb;
        
        Me = Me + w*detJac*me;
        Ke = Ke + w*detJac*ke;
        Ce = Ce + w*detJac*ce;
        Ve = Ve + w*detJac*ve;
    end
    % Store the locationa for assembly
	[i,j] = meshgrid(locu,locu);
	I(e,:) = i(:);
    J(e,:) = j(:);
    K(e,:) = locp;
    [x,y] = meshgrid(locu,locp);
    X(e,:) = x(:);
    Y(e,:) = y(:);
    
    M1(e,:) = reshape(Me,[numel(Me) 1]);
    K1(e,:) = reshape(Ke,[numel(Ke) 1]);
    C1(e,:) = reshape(Ce.',[numel(Ce) 1]);
    V1(e,:) = Ve;
end

% Preparing locations in 1D arrays for assembly
I = reshape(I,[numel(I) 1]);
J = reshape(J,[numel(J) 1]);
X = reshape(X,[numel(X) 1]);
Y = reshape(Y,[numel(Y) 1]);

% Preparing data in 1D arrays for assembly
M1 = reshape(M1,[numel(M1) 1]);
K1 = reshape(K1,[numel(K1) 1]);
C1 = reshape(C1,[numel(C1) 1]);

% Final assembly
Mmat = sparse(I,J,M1);
Kmat = sparse(I,J,K1);
Cmat = sparse(X,Y,C1);
Vmat = sparse(K,K,V1);
