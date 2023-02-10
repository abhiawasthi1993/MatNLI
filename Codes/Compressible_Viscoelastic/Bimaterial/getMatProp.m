function G = getMatProp(mesh)

X = mesh.P(:,1);
Y = mesh.P(:,2);
ly = mesh.ly;

for i = 1:numel(X)
    y = Y(i);
    if (y > (ly/2))
        G(i,1) = (10+1i*1)*1e-3;
    else
        G(i,1) = (20+1i*2)*1e-3;
    end
end