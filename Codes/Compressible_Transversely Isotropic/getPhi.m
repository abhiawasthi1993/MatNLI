function phi = getPhi(mesh)

X = mesh.P(:,1);
Y = mesh.P(:,2);

phi = zeros(numel(X),1);

for i = 1:numel(X)
    x = X(i);
    y = Y(i);

    t2 = 20 * exp(-(x-60)^2/(20^2) - (y-60)^2/(20^2));
    phi(i,1) = 0.01*(20+t2);
end
