function A = getFiberOrien(mesh, beta)

x1 = 1;     % Arbitrary x1 coordinate of x' vector
y2 = 1;     % Arbitrary y2 coordinate of y' vector

v1 = [x1 x1*tand(beta)];        % End coordinates of x' vector
v2 = [-v1(2)*(y2/x1) y2];       % End coordinates of y' vector (from dot(x',y')=0 => y1 = -x2*y2/x1)

a1 = v1./norm(v1);      % Unit vector originating from Origin (0,0)
a2 = v2./norm(v2);      % Unit vector originating from Origin (0,0)

a = [a1 a2];            % Local basis (e1' and e2')

A = ones(mesh.ne,4).*a; % Local basis for every element