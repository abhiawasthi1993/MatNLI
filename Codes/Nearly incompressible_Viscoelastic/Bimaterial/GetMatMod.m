function G = GetMatMod(mesh)

X = mesh.P(:,1);
Y = mesh.P(:,2);

for i = 1:numel(X)
    x = X(i);
    y = Y(i);
    
    if (y>60)
        G(i,1) = (10+1i*1)*1e-3;
    else
        G(i,1) = (20+1i*2)*1e-3;
    end
end