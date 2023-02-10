function stop = plotPar(E, optimValues, state)

if rem(optimValues.iteration,10) == 0
    load('data.mat','mesh');
    P = mesh.P;
    NCA = mesh.NCA;

    x = P(:,1);
    y = P(:,2);
    z = E;

    tri = [ NCA(1:end,[1 2 3]);  NCA(1:end,[1 3 4]) ];

    figure(2)
    trisurf(tri,x,y,z);
    xlabel('X-Axis')
    ylabel('Y-Axis')
    zlabel('Young''s Modulus')

    colormap('jet')
    shading interp
    colorbar
    view(0,90)
    drawnow
end

stop = false;