function ZS = FitMLSSurface(X,Y,Z,deltaX,deltaY)

noPts = size(X,1);
node = [X Y];
numnode = size(X,1);

shape = 'circle' ;         % shape of domain of influence
dmax  = 3 ;                % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function

delta  = max(deltaX,deltaY);
di     = ones(1,numnode)*dmax*delta ;

parfor i=1:noPts
    pt = [X(i) Y(i)];
    [index] = define_support(node,pt,di);
    [phi,dphidx,dphidy] = MLS_ShapeFunction(pt,index,node,di,form);
    ZS(i) = phi*Z(index);
end
end