function Gmat = getAdjointPar(Gvec,mesh,Nu,U,Pr)

P = mesh.P;
NCA = mesh.NCA;
nu = size(P,1);
ne = size(NCA,1);
np = ne;

ngp = 2;                % Number of Gauss points
[GP,W] = Gauss(ngp);    % Gauss point locations and weights

% Gauss point locations and weights
parfor e = 1:ne
    a = NCA(e,:);
    xe = P(a,1); ye = P(a,2);
    locu = zeros(1,8);
    locu(1:2:8) = 2*a-1;
    locu(2:2:8) = 2*a;
    locp = a;

    % Get nodal quantities for element 'e'
    Ge = Gvec(locp,1);         % Complex shear modulus
    Ue = U(locu,1);         % Displacement field
    Pre = Pr(e,1);          % Pressure field
    
    dKdGe = zeros(8,4);
    dVdGe = zeros(1,4);

    % Computation at Gauss point level
    for g = 1:size(GP,1)
        r = GP(g,1);
        s = GP(g,2);
        w = W(g);
        
        % Get shape functions and their derivatives
        [N,DNr,DNs] = getQ4shp(r,s);
        Jac = [DNr; DNs]*[xe ye];
        detJac = det(Jac);
        DN = Jac\[DNr;DNs];
        
        Nd = zeros(2,8);
        Nd(1,1:2:8) = N;
        Nd(2,2:2:8) = N;

        % Compute strain-displacement matrix (B)
        B = zeros(3,8);
        B(1,1:2:8) = DN(1,:);
        B(2,2:2:8) = DN(2,:);
        B(3,1:2:8) = DN(2,:);
        B(3,2:2:8) = DN(1,:);
        
        % Interpolate nodal quantities to the Gauss points
        E = 2*(N*Ge)*(1+Nu);        % Young's modulus
        mu = (N*Ge);                % Shear modulus
        coef = 1/(3*(1-2*Nu));
        Kb = E*coef;                % Bulk modulus
		
		% Deviatoric part of constitutive matrix
        I0 = 0.5*[2 0 0;0 2 0; 0 0 1];
        m = [1 1 0]';
        Dd = 2*(I0 - 1/3*m*m');
		
        % Compute matrices
		ge = zeros(8,4);
        ge = B'*Dd*(B*Ue)*N;
		he = zeros(1,4);
        he = 1*Pre*(1/Kb^2)*coef*N*2*(1+Nu);

        dKdGe = dKdGe + w*detJac*ge;
		dVdGe = dVdGe + w*detJac*he;
    end
    % Store the locations for assembly
    [i,j] = meshgrid(locu,locp);
    I(e,:) = i(:);
    J(e,:) = j(:);
    [x,y] = meshgrid(e,locp);
    X(e,:) = x(:);
    Y(e,:) = y(:);
    dKdG(e,:) = reshape(dKdGe.',[numel(dKdGe) 1]);
	dVdG(e,:) = reshape(dVdGe.',[numel(dVdGe) 1]);
end

% Preparing locations in 1D arrays for assembly
I = reshape(I, [numel(I) 1]);
J = reshape(J, [numel(J) 1]);
X = reshape(X,[numel(X) 1]);
Y = reshape(Y,[numel(Y) 1]);

% Preparing data in 1D arrays for assembly
dKdG = reshape(dKdG, [numel(dKdG) 1]);
dVdG = reshape(dVdG, [numel(dVdG) 1]);

% Final assembly
Gmat1 = sparse(I,J,dKdG);
Gmat2 = sparse(X,Y,dVdG);

Gmat = [Gmat1; Gmat2];