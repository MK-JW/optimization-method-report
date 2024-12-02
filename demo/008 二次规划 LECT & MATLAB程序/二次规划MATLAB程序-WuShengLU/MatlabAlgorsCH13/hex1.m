% Program: hex1.m
% Description: This function can be used to 
% test function kelley_ie.
 function h = hex1(x)
 xx = x(1:2);
 H = [4 1; 1 1];
 p = [-1 2]';
 h1 = [H*xx+p; -1];
 h2 = [1 -1 0]';
 h3 = [1 1 0]';
 h4 = [-1 0 0]';
 h5 = [0 -1 0]';
 h = [h1 h2 h3 h4 h5];