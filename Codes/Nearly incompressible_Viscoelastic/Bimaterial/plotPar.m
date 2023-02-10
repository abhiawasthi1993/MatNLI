function stop = plotPar(E, optimValues, state)

if rem(optimValues.iteration,10) == 0
    load('data.mat','mesh');

    P = mesh.P;
    nn = mesh.nn;
    NCA = mesh.NCA;

    x=P(:,1);
    y=P(:,2);
    z=E(1:nn,1) + 1i*E(nn+[1:nn],1);

    tri = [ NCA(1:end,[1 2 3]);  NCA(1:end,[1 3 4]) ];

    figure(2)
    subplot(1,2,1)
    trisurf(tri,x,y,real(z));
    xlabel('X-Axis')
    ylabel('Y-Axis')
    zlabel('Young''s Modulus')
    axis tight equal
    colormap('jet')
    shading interp
    colorbar
    view(0,90)
    drawnow

    figure(2)
    subplot(1,2,2)
    trisurf(tri,x,y,imag(z));
    xlabel('X-Axis')
    ylabel('Y-Axis')
    zlabel('Young''s Modulus')
    axis tight equal
    colormap('jet')
    shading interp
    colorbar
    view(0,90)
    drawnow
end

stop = false;