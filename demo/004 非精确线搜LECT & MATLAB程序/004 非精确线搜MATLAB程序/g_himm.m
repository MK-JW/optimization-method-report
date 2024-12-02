% g_himm.m
function z = g_himm(x)
x1 = x(1);
x2 = x(2);
w1 = (x1^2 + x2 - 11);
w2 = (x1 + x2^2 - 7);
z1 = 4*w1*x1 + 2*w2;
z2 = 2*w1 + 4*w2*x2;
z = [z1 z2]';