function [GP,W] = Gauss(ngp)

if ngp == 2
    q = sqrt(1/3);
    GP = [-q -q;
           q -q;
          -q  q;
           q  q];
    W = [1 1 1 1];
end

if ngp == 3
    q = sqrt(3/5);
    w1 = 5/9;
    w2 = 8/9;
    GP = [ -q -q;
            0 -q
            q -q;
           -q  0;
            0  0;
            q  0;
           -q  q;
            0  q;
            q  q ];

    W = [w1*w1;
         w2*w1;
         w1*w1;
         w1*w2;
         w2*w2;
         w1*w2;
         w1*w1;
         w2*w1;
         w1*w1];
end
