% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                            Nguyen Vinh Phu
%                        LTDS, ENISE, Juillet 2006
% *************************************************************************
function [index] = define_support(node,x,di)
% find nodes in neighbouring of point x
% Inpputs:
%  node : numnode x 2, nodal coordinates
%  x    : 1 x 2, coordinate of point 
%  di   : 1 x numnode, size of support of nodes

numnode = size(node,1) ;
dif = node - [ones(numnode,1)*x(1,1) ones(numnode,1)*x(1,2)];
for i = 1 : numnode
    r(i) = norm(dif(i,:));
end

index = find(r - di <= 0.00001);
length(index);