function plotPress(P,NCA,Pr,fno,titlefig)

x=P(:,1);
y=P(:,2);

z = 0*x;
c = 0*x;

z(NCA(:,1),1) = z(NCA(:,1),1) + Pr;
z(NCA(:,2),1) = z(NCA(:,2),1) + Pr;
z(NCA(:,3),1) = z(NCA(:,3),1) + Pr;
z(NCA(:,4),1) = z(NCA(:,4),1) + Pr;

c(NCA(:,1),1) = c(NCA(:,1),1) + 1;
c(NCA(:,2),1) = c(NCA(:,2),1) + 1;
c(NCA(:,3),1) = c(NCA(:,3),1) + 1;
c(NCA(:,4),1) = c(NCA(:,4),1) + 1;

z = z./c;

tri1 = NCA(1:end,[1 2 3]);
tri2 = NCA(1:end,[1 3 4]);
tri =[tri1; tri2];

figure(fno)
hold on
trisurf(tri,x,y,z);

xlabel('X-Axis')
ylabel('Y-Axis')
zlabel(titlefig)
colormap('jet')
shading interp
colorbar
title(titlefig)
view(0,90)
set(gcf,'color','w');
set(gca,'FontSize',30)
axis equal tight
