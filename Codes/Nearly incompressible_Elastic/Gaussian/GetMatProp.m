function E = GetMatProp(mesh)

X = mesh.P(:,1);
Y = mesh.P(:,2);

for i = 1:numel(X)
    x = X(i);
    y = Y(i);
    t2 = 5*exp(-(x-60)^2/(20^2) - (y-60)^2/(20^2));
    E(i,1) = 0.001*(10+t2);
end