function chi = getZeta(mesh)

X = mesh.P(:,1);
Y = mesh.P(:,2);

chi = zeros(numel(X),1);

for i = 1:numel(X)
    x = X(i);
    y = Y(i);
    t2 = 20 * exp(-(x-60)^2/(20^2) - (y-60)^2/(20^2));
    chi(i,1) = 0.1*(5 + t2);
end
