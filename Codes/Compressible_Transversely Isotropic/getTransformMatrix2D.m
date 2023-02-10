function T = getTransformMatrix2D(a)

a1 = a(1,1:2); a2 = a(1,3:4);   % Local basis (e1' and e2')

e1 = [1 0];    e2 = [0 1];      % Global basis (e1 and e2)

% Direction cosines
t11 = dot(a1,e1);  t12 = dot(a1,e2);
t21 = dot(a2,e1);  t22 = dot(a2,e2);

% Transformation matrix
T = [ t11*t11     t12*t12        t11*t12;
      t21*t21     t22*t22        t21*t22;
     2*t11*t21   2*t12*t22  (t11*t22 + t12*t21)];