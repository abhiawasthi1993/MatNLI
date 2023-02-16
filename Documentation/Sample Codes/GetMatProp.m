function G = GetMatProp(mesh)

X = mesh.P(:,1);
Y = mesh.P(:,2);

xc = 60;            % x-coordinate of Gaussian center
yc = 60;            % y-coordinate of Gaussian center

for i = 1:numel(X)
    x = X(i);
    y = Y(i);
    r2 = (x-xc)^2 + (y-yc)^2;
    R = 10 + 5*exp(-r2/400);
    C = 1 + 0.5*exp(-r2/400);
    G(i,1) = R + 1i*C;
end

G = G*1e-3;     % Convert into MPa
