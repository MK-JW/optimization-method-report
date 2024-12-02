% Program: fex.m
% Description: This function can be used to 
% test function ellipsoid_ie.
function f = fex(x)
 H = [4 1; 1 1];
 p = [-1 2]';
 f = 0.5*x'*H*x + x'*p;