close all;
clear all;
clc;
figure
hold on

kk=1:1:50;

x0=[1,3.5]';
[x,mu,lambda,output,fx_rec]=multphr_figure('f1','h1','g1','df1','dh1','dg1',x0);
n=length(fx_rec);
plot(kk(1:n),fx_rec(1:n),'color','b','Marker','.');
plot(kk(n),fx_rec(n),'color','b','Marker','p')

x0=[3,3]';
[x,mu,lambda,output,fx_rec]=multphr_figure('f1','h1','g1','df1','dh1','dg1',x0);
n=length(fx_rec);
plot(kk(1:n),fx_rec(1:n),'color','m','Marker','.');
plot(kk(n),fx_rec(n),'color','m','Marker','p')


x0=[1,1]';
[x,mu,lambda,output,fx_rec]=multphr_figure('f1','h1','g1','df1','dh1','dg1',x0);
n=length(fx_rec);
plot(kk(1:n),fx_rec(1:n),'color','r','Marker','.');
plot(kk(n),fx_rec(n),'color','r','Marker','p')


x0=[2,4]';
[x,mu,lambda,output,fx_rec]=multphr_figure('f1','h1','g1','df1','dh1','dg1',x0);
n=length(fx_rec);
plot(kk(1:n),fx_rec(1:n),'color','c','Marker','.');
plot(kk(n),fx_rec(n),'color','c','Marker','p')