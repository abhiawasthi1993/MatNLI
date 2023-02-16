function [w,dwdx,dwdy] = circle_spline(x,xI,d,form)
% Compute cubic and quartic spline function
% Inputs:
% x (1x2)  : coordinate of point at which w is to be evaluated
% xI (1x2) : coord of node I
% d        : size of the support

r = sqrt( (x(1,1) - xI(1,1)) .* (x(1,1) - xI(1,1)) + (x(1,2) - xI(1,2)).*(x(1,2) - xI(1,2)) )/d ;

switch form
  case 'cubic_spline' 
     [w,dwdr] = cubic_spline(r);
  case 'quartic_spline'
     [w,dwdr] = quartic_spline(r);
  otherwise 
     error('Grr. Unknown functional form');
end

if (r ~= 0)
    drdx = (x(1,1) - xI(1,1))/(r*d*d) ;
    drdy = (x(1,2) - xI(1,2))/(r*d*d) ;
else
    drdx = 0 ;
    drdy = 0 ;
end

dwdx = dwdr * drdx ;
dwdy = dwdr * drdy ;

