function Example_13_5_Jacob_form
% Find the shortest distance between the two ellipses using Algorithm13.6 |
% Algorithm 13.6 Kelley s cutting-plane algorithm for CP problems         |
% with inequality constraints                                             |
%==========================================================================
% min 0.5*x'*H*x                   % distance between ellipse 1 and 2
% s.t. 0.5*x'*Q1*x+x'*q1+r1 >= 0   % ellipse 1
%      0.5*x'*Q2*x+x'*q2+r2 >= 0   % ellipse 2
% H=[1,0,-1,0;0,1,0,-1;-1,0,1,0;0,-1,0,1];
% Q1=-[0.5,0,0,0;0,2,0,0;0,0,0,0;0,0,0,0]; q1=[0.5;0;0;0]; r1=3/4;
% Q2=-[0,0,0,0;0,0,0,0;0,0,5/4,3/4;0,0,3/4,5/4]; q2=[0;0;11/2;13/2]; 
% r2=-35/2;
%========================================================================
close all;clear all;clc;
format long
global x_rec
%initial point X0=[x0,L0]'
%x0=[1.5,0.5,2.5,4.0]';
x0=[1.5,0.5,2.5,3.0]';
L0=4;
X0=[x0;L0];
%===========n=5,p=3 so we consider another 2 initial points============
%           then we have 3 initial points and 3*3=9 constraints
% because for LP, we have optimal solution while p>=n  
X1=[0,0,1,5,16]';
X2=[2,0,2,4,9]';
%===========keep Xi feasible for L-f(x)>=0 ==================
c0=cex_13_5(X0);
while c0(1)<0
    X0(5)=c0(1)+0.5;
    c0=cex_13_5(X0);
end
%
c1=cex_13_5(X1);
while c1(1)<0
    X1(5)=c1(1)+0.5;
    c1=cex_13_5(X1);
end
%
c2=cex_13_5(X2);
while c2(1)<0
    X2(5)=c2(1)+0.5;
    c2=cex_13_5(X2);
end
%======================================================
X0=[X0,X1,X2];

epsi = 1e-7;

%============================================================
[xs,fs,k]= kelley_ie_d('cex_13_5','dex_13_5',X0,epsi)

r=xs(1:2)
s=xs(3:4)
dist=norm(r-s,2)
%=============Figure===============================================
figure(1)
ezplot('-0.25*tx.^2-ty.^2+0.5*tx+0.75')
hold on
ezplot('-(5*tx.^2+5*ty.^2+6*tx.*ty)/8+5.5*tx+6.5*ty-17.5')
axis([-2 4.5 -2.5 6])

xtext=1;ytext=4.0;
text(xtext,ytext,'S','FontSize',16)
xtext=0;ytext=0;
text(xtext,ytext,'R','FontSize',16)

plot(x0(1),x0(2),'or')
plot(x0(3),x0(4),'or')
plot(X1(1),X1(2),'or')
plot(X1(3),X1(4),'or')
plot(X2(1),X2(2),'or')
plot(X2(3),X2(4),'or')

plot(xs(1),xs(2),'ok')
plot(xs(3),xs(4),'ok')

%======path========================

hold on
%plot(x_rec(:,1),x_rec(:,2),'g')
plot(x_rec(:,1),x_rec(:,2),'.g')
%plot(x_rec(:,3),x_rec(:,4),'b')
plot(x_rec(:,3),x_rec(:,4),'.b')

plot([r(1),s(1)],[r(2),s(2)],'--k')

%==========================================
plot(r(1),r(2),'*r')
plot(s(1),s(2),'*r')
hold on
xtext=r(1)+0.20;ytext=r(2)+0.2;
text(xtext,ytext,'r^*','Color','m','FontSize',12)
xtext=s(1)+0.0;ytext=s(2)-0.4;
text(xtext,ytext,'s^*','Color','m','FontSize',12)

str=['Distance between R and S =',num2str(dist)]
title(str)
xlabel('x_1 / x_3');ylabel('x_2 / x_4');

%===========figure Xi=========================
figure(2)
[np,mp]=size(x_rec)
subplot(3,2,1)
plot([0:mp-1],x_rec(1,:))
hold on
xlabel('Num of Iter');ylabel('x_1');
subplot(3,2,2)
plot([0:mp-1],x_rec(2,:))
hold on
xlabel('Num of Iter');ylabel('x_2');
subplot(3,2,3)
plot([0:mp-1],x_rec(3,:))
hold on
xlabel('Num of Iter');ylabel('x_3');
subplot(3,2,4)
plot([0:mp-1],x_rec(4,:))
hold on
xlabel('Num of Iter');ylabel('x_4');
fs=sqrt(((x_rec(1,:)-x_rec(3,:)).^2+(x_rec(2,:)-x_rec(4,:)).^2));
subplot(3,2,[5,6])
plot([0:mp-1],fs)
hold on
xlabel('Num of Iter');ylabel('dist');
end


%================cex_13_5.m============================================
% Program: cex_13_5.m
% Description: This function can be used to 
% test function kelley_ie.
% x=[x_points;L],example: x=X0=[x0;L0]
 function c = cex_13_5(x)
 format long
 nx=length(x);
 L = x(end);
 xx = x(1:4);
 x1=xx(1:2);
 x2=xx(3:4);
 H=[1,0,-1,0;0,1,0,-1;-1,0,1,0;0,-1,0,1];
 Q1=-[0.5,0;0,2];
 Q2=-[5/4,3/4;3/4,5/4];
 q1=[0.5;0];q2=[11/2;13/2];
 r1=3/4;r2=-35/2;
 c1 = L - (0.5*xx'*H*xx );
 c2 = 0.5*x1'*Q1*x1 +x1'*q1+r1;
 c3 = 0.5*x2'*Q2*x2 +x2'*q2+r2;
 c = [c1;c2;c3];
 end

 %================dex_13_5.m============================================
 % Program: dex_13_5.m
% Description: This function can be used to 
% test function kelley_ie.  Jacobian matrix
function d = dex_13_5(x)
 format long
 nx=length(x);
 L = x(nx);
 xx = x(1:4);
 x1=xx(1:2);
 x2=xx(3:4);
 H=[1,0,-1,0;0,1,0,-1;-1,0,1,0;0,-1,0,1];
 Q1=-[0.5,0;0,2];
 Q2=-[5/4,3/4;3/4,5/4];
 q1=[0.5;0];q2=[11/2;13/2];


 d1 = [-xx'*H, 1];
 d2 = [x1'*Q1+q1',0,0,0];
 d3 = [0,0,x2'*Q2+q2',0];
 
 d = [d1;d2;d3];
end

%==============Kelley=s cutting-plane algorithm===========================
% Program: kelley_ie.m
% Title: Kelley's cutting-plane algorithm for the CP 
% problem in Eq, (13.70).
% Description: Implements Algorithm 13.6 for the CP 
% problem in Eq. (13.70).
% Theory: See Practical Optimization Sec. 13.5.
% Input:     
%         cname -- name of MATLAB function that evaluates
%                  functions cj(x)for j = 1, 2, ..., q in 
%                  Eq.(13.70b).
%         hname -- name of MATLAB function that evaluates
%                  gradients of -cj(x) for j = 1,2,...,q.
%            X0 -- matrix whose columns are a set of 
%                  feasible points
%          epsi -- convergence tolerance
% Output:   
%            xs -- solution vector
%            fs -- value of objective function at xs
%             k -- number of iterations at convergence
% Example:
% Apply Algorithm 13.6 to solve the QP problem
% min 0.5*x'*H*x                   % distance between ellipse 1 and 2
% s.t. 0.5*x'*Q1*x+x'*q1+r1 >= 0   % ellipse 1
%      0.5*x'*Q2*x+x'*q2+r2 >= 0   % ellipse 2
% H=[1,0,-1,0;0,1,0,-1;-1,0,1,0;0,-1,0,1];  
% Q1=-[0.5,0,0,0;0,2,0,0;0,0,0,0;0,0,0,0]; q1=[0.5;0;0;0]; r1=3/4;
% Q2=-[0,0,0,0;0,0,0,0;0,0,5/4,3/4;0,0,3/4,5/4]; q2=[0;0;11/2;13/2]; 
% r2=-35/2;
%
% Solution:
% First, convert the problem at hand to the form of
% 
% Function  kelley_ie requires two MATLAB functions cex1.m 
% and hex1.m as follows:
%
% =======================
% These m-files can be found in the same folder as m-file kelley_ie.m.
% The solution can be obtained by executing the following commands:
% X0 = [1 1 2 2; 1 2 1 2; 5 10 11 18]
% epsi = 1e-7
% [xs,fs,k]= kelley_ie_d('cex1','dex1',X0,epsi)
% =================================================
function [xs,fs,k] = kelley_ie_d(cname,dname,X0,epsi)
format long
global x_rec
NMAX=100;
disp(' ')
disp('Program kelley_ie.m')
[n,K] = size(X0);
Ak = [];
bk = [];
c = [zeros(n-1,1);1];
x0 = X0(:,1);
k_rec=1;
x_rec=[];
x_rec(:,k_rec)=x0(1:4);
for i = 1:K
    xi = X0(:,i);
    ci = feval(cname,xi);
    di = feval(dname,xi);
    Ai = di;
    Ak = [Ak; Ai];
    bk = [bk; Ai*xi-ci];
end

k = 1;
%[xw,fsw,kw] = lp_ad(Ak,bk,c,x0)
[xw,fsw]=linprog(c,-Ak,-bk);
k_rec=k_rec+1;
x_rec(:,k_rec)=xw(1:4);
cw = feval(cname,xw);
[yw,iw] = min(cw);
KK=0;
while yw < -epsi && KK < NMAX
    KK=KK+1;
    cjs = cw(iw);
    d1 = feval(dname,xw);
    djs = d1(iw,:);
    Ak = [Ak; djs];
    bk = [bk; djs*xw-cjs];
    %[xw,fsw,kw] = lp_ad(Ak,bk,c,x0);
    [xw,fsw]=linprog(c,-Ak,-bk);
    k_rec=k_rec+1;
    x_rec(:,k_rec)=xw(1:4);
    cw = feval(cname,xw);
    [yw,iw] = min(cw);
    k = k + 1;
end
xs = xw;
fs = xs(end); 
end