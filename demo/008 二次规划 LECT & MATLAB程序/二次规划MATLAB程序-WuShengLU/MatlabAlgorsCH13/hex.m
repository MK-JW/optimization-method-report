% Program: hex.m
% Description: This function can be used to 
% test function ellipsoid_ie.
 function h = hex(x)
 h1 = [1 -1]';
 h2 = [1 1]';
 h3 = [-1 0]';
 h4 = [0 -1]';
 h = [h1 h2 h3 h4];