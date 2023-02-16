function E = getMatProp(mesh)

X = mesh.P(:,1);
Y = mesh.P(:,2);

for i = 1:numel(X)
    x = X(i);
    y = Y(i);
    t1 = 15 * exp(-(x-10)^2/6.6 - (y-5)^2/6.6);
    t2 = 30 * exp(-(x-10)^2/6.6 - (y-15)^2/6.6);
    E(i,1) = 0.001*(10 + t1 + t2);
end
