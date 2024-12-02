% Program: cex1.m
% Description: This function can be used to 
% test function kelley_ie.
 function c = cex1(x)
 x1 = x(1);
 x2 = x(2);
 L = x(3);
 xx = x(1:2);
 H = [4 1; 1 1];
 p = [-1 2]';
 c1 = L - (0.5*xx'*H*xx +xx'*p);
 c2 = -x1 + x2 + 2;
 c3 = -x1 - x2 + 6;
 c4 = x1;
 c5 = x2;
 c = [c1 c2 c3 c4 c5]';