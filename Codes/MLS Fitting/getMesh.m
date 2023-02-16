function mesh = getMesh(lx,ly,nx,ny)

x = linspace(0,lx,nx);
y = linspace(0,ly,ny);

z = 0;
for j = 1:numel(y)
    for i = 1:numel(x)
        z = z+1;
        P(z,1:2) = [x(i) y(j)];
    end
end

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

ln = 1:nx:(nx*ny);
rn = nx:nx:(nx*ny);
bn = 1:nx;
tn = (nx*ny-nx+1):(nx*ny);

BC = unique([ln rn bn tn]);
hDBC = [2*ln-1 2*tn];
DBC = [2*bn];

mesh.nn = size(P,1);
mesh.ne = size(NCA,1);
mesh.P = P;
mesh.NCA = NCA;
mesh.hDBC = hDBC;
mesh.DBC = DBC;
mesh.BC = BC;
mesh.LB = ln;
mesh.RB = rn;
mesh.TB = tn;
mesh.BB = bn;
mesh.lx = lx;
mesh.ly = ly;
mesh.nx = nx;
mesh.ny = ny;
