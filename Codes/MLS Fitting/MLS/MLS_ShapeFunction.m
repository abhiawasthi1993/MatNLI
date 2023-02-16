function [phi,dphidx,dphidy] = MLS_ShapeFunction(pt,index,node,di,form)

% Compute the MLS shape function at point pt for all nodes within the
% support of this point pt.
% Basis used is linear basis pT = [1 x y]

% --------------------------------------
%      compute the moment matrix A
% --------------------------------------

A    = zeros(3,3) ;
dAdx = zeros(3,3) ;
dAdy = zeros(3,3) ;
for m = 1 : size(index,2)
    xi = [node(index(m),1) node(index(m),2)] ;
    [wi,dwidx,dwidy] = circle_spline(pt,xi,di(index(m)),form);
    pTp = [1 xi(1,1) xi(1,2)]'*[1 xi(1,1) xi(1,2)] ; % Linear basis, pT = [1 x y]
    A    = A    + wi*pTp ;
    dAdx = dAdx + dwidx*pTp ;
    dAdy = dAdy + dwidy*pTp ;
    % store weight function and its derivative at node I for later use
    w(m)    = wi ;
    dwdx(m) = dwidx ;
    dwdy(m) = dwidy ;
end

clear wi; clear dwidx; clear dwidy ; clear xi;

p  = [1; pt(1,1); pt(1,2)];

% --------------------------------------
%         compute  matrix c(x)
% --------------------------------------

% A(x)c(x)   = p(x)
% A(x)c,k(x) = b,k(x)
% Backward substitutions, two times for c(x), two times for c,k(x) k
% =1,2

% Using LU factorization for A
[L,U,PERM] = lu(A) ;

for i = 1 : 3
    if i == 1         % backward substitution for c(x)
        C = PERM*p;
    elseif i == 2     % backward substitution for c,x(x)
        C = PERM*([0 1 0]' - dAdx*c(1:3,1));
    elseif i == 3     % backward substitution for c,y(x)
        C = PERM*([0 0 1]' - dAdy*c(1:3,1));
    end

    D1 = C(1);
    D2 = C(2) - L(2,1)*D1;
    D3 = C(3) - L(3,1)*D1 - L(3,2)*D2 ;

    c(3,i) = D3/U(3,3) ;
    c(2,i) = (D2 - U(2,3)*c(3,i))/(U(2,2));
    c(1,i) = (D1 - U(1,2)*c(2,i) - U(1,3)*c(3,i))/(U(1,1));
end

for m = 1 : size(index,2)
    xi = [node(index(m),1) node(index(m),2)] ;
    piT = [1 xi(1,1) xi(1,2)]';
    phi(m) = c(:,1)'* piT*w(m) ;
    dphidx(m) = c(:,2)'*piT*w(m) + c(:,1)'*piT*dwdx(m) ;
    dphidy(m) = c(:,3)'*piT*w(m) + c(:,1)'*piT*dwdy(m);
end
