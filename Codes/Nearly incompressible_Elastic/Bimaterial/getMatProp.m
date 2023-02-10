function E = getMatProp(mesh)

X = mesh.P(:,1);
Y = mesh.P(:,2);
ly = mesh.ly;

for i = 1:numel(X)
    y = Y(i);
    if (y > (ly/2))
        E(i,1) = (10)*1e-3;
    else
        E(i,1) = (20)*1e-3;
    end
end
