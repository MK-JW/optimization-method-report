function test0
x=[4;-1];
hval = hfun(x)
for i=1:4
    hval =hfun(x,i)
end
end

function hval = hfun(x,i)
% c1 =  2*x1-x2- 6 >= 0
% c2 = -2*x1-x2+10 >= 0
% c3 = -3*x1+x2+15 >= 0
% c4 =  3*x1+x2- 9 >= 0
% h =[-▽c1,-▽c2,-▽c3,-▽c4]
hval=-[  2, -2, -3,  3;
        -1, -1,  1,  1];
if nargin < 2
    hval = hval;
else
    hval = hval(:,i);
end
end

