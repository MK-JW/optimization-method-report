% Program: gex.m
% Description: This function can be used to 
% test function ellipsoid_ie.
function g = gex(x)
 H = [4 1; 1 1];
 p = [-1 2]';
 g = H*x + p;