% Program: cex.m
% Description: This function can be used to 
% test function ellipsoid_ie.
 function c = cex(x)
 x1 = x(1);
 x2 = x(2);
 c1 = -x1 + x2 + 2;
 c2 = -x1 - x2 + 6;
 c3 = x1;
 c4 = x2;
 c = [c1 c2 c3 c4]';