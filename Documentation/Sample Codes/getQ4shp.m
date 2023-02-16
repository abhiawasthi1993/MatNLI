function [N,DNr,DNs] = getQ4shp(r,s)

N1 = 0.25*(1-r)*(1-s);
N2 = 0.25*(1+r)*(1-s);
N3 = 0.25*(1+r)*(1+s);
N4 = 0.25*(1-r)*(1+s);

N1r = 0.25*(-1)*(1-s);
N2r = 0.25*(+1)*(1-s);
N3r = 0.25*(+1)*(1+s);
N4r = 0.25*(-1)*(1+s);

N1s = 0.25*(1-r)*(-1);
N2s = 0.25*(1+r)*(-1);
N3s = 0.25*(1+r)*(+1);
N4s = 0.25*(1-r)*(+1);

N   = [N1 N2 N3 N4];        % Shape function
DNr = [N1r N2r N3r N4r];    % Derivative of shape function w.r.t. 'r'
DNs = [N1s N2s N3s N4s];    % Derivative of shape function w.r.t. 's'




