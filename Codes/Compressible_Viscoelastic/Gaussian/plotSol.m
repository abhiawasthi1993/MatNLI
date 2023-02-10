function plotSol(P,NCA,U,fno,titlefig)

x=P(:,1);
y=P(:,2);
z=U;

tri = [ NCA(1:end,[1 2 3]);  NCA(1:end,[1 3 4]) ];
    
figure(fno)
trisurf(tri,x,y,z);
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel(titlefig)

colormap('jet')
shading interp
colorbar
title(titlefig)
view(0,90)
drawnow
set(gcf,'color','w');
set(gca,'FontSize',30)
axis equal tight
