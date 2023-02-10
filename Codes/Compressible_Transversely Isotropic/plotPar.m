function stop = plotPar(E1, optimValues, state)

load('ScData.mat','D','Ct');
itr = optimValues.iteration;

if rem(itr,1000)==0
    E = D*E1 + Ct;
    load('data.mat','mesh','Gvec');
    P = mesh.P;
    NCA = mesh.NCA;
    nn = mesh.nn;

    x = P(:,1);
    y = P(:,2);

    tri = [ NCA(1:end,[1 2 3]);  NCA(1:end,[1 3 4]) ];

    figure(2)
    subplot(2,2,1)
    trisurf(tri,x,y,E(1:nn,1));
    xlabel('X-Axis')
    ylabel('Y-Axis')
    title('\mu')
    colormap('jet')
    shading interp
    colorbar
    view(0,90)
    axis tight equal
    drawnow

    figure(2)
    subplot(2,2,2)
    trisurf(tri,x,y,E(nn+[1:nn],1));
    xlabel('X-Axis')
    ylabel('Y-Axis')
    title('\zeta')
    colormap('jet')
    shading interp
    colorbar
    view(0,90)
    axis tight equal
    drawnow

    figure(2)
    subplot(2,2,3)
    trisurf(tri,x,y,Gvec(1:nn,1));
    xlabel('X-Axis')
    ylabel('Y-Axis')
    title('\mu - Truth')
    colormap('jet')
    shading interp
    colorbar
    view(0,90)
    axis tight equal
    drawnow

    figure(2)
    subplot(2,2,4)
    trisurf(tri,x,y,E(2*nn+[1:nn],1));
    xlabel('X-Axis')
    ylabel('Y-Axis')
    title('\phi')
    colormap('jet')
    shading interp
    colorbar
    view(0,90)
    axis tight equal
    drawnow
end

stop = false;