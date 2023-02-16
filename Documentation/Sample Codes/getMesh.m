function mesh = getMesh(lx,ly,nx,ny)

x = linspace(0,lx,nx);
y = linspace(0,ly,ny);

z = 0;

% Create nodes and their coordinates (P)
for j=1:numel(y)
    for i=1:numel(x)
        z=z+1;
        P(z,1:2) = [x(i) y(j)];
    end
end

% Form element connectivity (NCA)
e = 0;
for j = 1:(numel(y)-1)
    for i = 1:(numel(x)-1)
        e = e+1;
        a1 = nx*(j-1) + i;
        a2 = a1 + 1;
        a3 = a2 + nx;
        a4 = a1 + nx;
        NCA(e,1:4) = [a1 a2 a3 a4];
    end
end

% Identify boundary nodes
ln = 1:nx:(nx*ny);          % Nodes on left boundary
rn = nx:nx:(nx*ny);         % Nodes on right boundary
bn = 1:nx;                  % Nodes on bottom boundary
tn = (nx*ny-nx+1):(nx*ny);  % Nodes on top boundary

BC = unique([ln rn bn tn]); % All boundary nodes
DBC = 2*tn;                 % y-DoFs corresponding to top boundary
NBC = 2*tn-1;               % x-DoFs corresponding to top boundary

% Save all data as a MATLAB structure array
mesh.nn = size(P,1);        % Number of nodes
mesh.ne = size(NCA,1);      % Number of elements
mesh.P = P;                 % Nodal coordinates
mesh.NCA = NCA;             % Nodal connectivity array
mesh.DBC = DBC;
mesh.NBC = NBC;
mesh.BC = BC;
mesh.LB = ln;
mesh.RB = rn;
mesh.TB = tn;
mesh.BB = bn;
