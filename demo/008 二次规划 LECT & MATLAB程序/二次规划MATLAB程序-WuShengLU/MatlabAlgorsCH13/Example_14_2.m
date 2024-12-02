function Example_14_2
close;clear;clc;
global y_rec
format long
[AM,X0,y0,S0,C,b,gama,epsl] = data_ex14_2;
Ag=AM;gam=gama;epsi=epsl;
[X,y,S,k] = sdp_pc(X0,y0,S0,Ag,b,C,gam,epsi)
r=[y(1);y(2)]
s=[y(3);y(4)]
dist=norm(r-s,2)

figure(1)
ezplot('-0.25*tx.^2-ty.^2+0.5*tx+0.75')
hold on
ezplot('-(5*tx.^2+5*ty.^2+6*tx.*ty)/8+5.5*tx+6.5*ty-17.5')
axis([-1.5 4 -1.5 6])

xtext=1;ytext=4.5;
text(xtext,ytext,'S','FontSize',16)
xtext=0;ytext=0;
text(xtext,ytext,'R','FontSize',16)

plot(y_rec(:,1),y_rec(:,2),'r')
plot(y_rec(:,1),y_rec(:,2),'.r')
plot(y_rec(:,3),y_rec(:,4),'k')
plot(y_rec(:,3),y_rec(:,4),'.k')

plot([y(1),y(3)],[y(2),y(4)],'g')

plot(y0(1),y0(2),'or')
plot(y0(3),y0(4),'ok')
plot(y(1),y(2),'*r')
plot(y(3),y(4),'*r')
hold on
xtext=y(1)+0.08;ytext=y(2)+0.2;
text(xtext,ytext,'r^*','Color','red','FontSize',12)
xtext=y(3)+0.05;ytext=y(4)-0.2;
text(xtext,ytext,'s^*','Color','red','FontSize',12)

str=['Distance between R and S =',num2str(dist)]
title(str)
xlabel('x_1 / x_3');ylabel('x_2 / x_4');
end

function [AM,X0,y0,S0,C,b,gama,epsl] = data_ex14_2
% Program: data_ex14_2.m
% Title: Data matrices of Example 14.2.
% Description: Loads the data matrices of Example 14.3 into the 
% MATLAB environment and generates matrices Ai,C and b for
% A strictly feasible set {X0,y0,S0}
% functions Predictor-Corrector Method.
% Input: None.     
% Output:    
%   Matrices AM(Ai),C,b and initial set {X0,y0,S0}
% Example:
% Execute
%   [AM,X0,y0,S0,C,b,gama,epsl] = data_ex14_2
% =================================================
H=[1,0,-1,0;0,1,0,-1;-1,0,1,0;0,-1,0,1];
Q1=[0.25,0;0,1];Q2=[5,3;3,5];
q1=[-0.5;0];q2=[-44;-52];
r1=-0.75;r2=140;
HH=[1,0,-1,0;0,-1,0,1;0,0,0,0;0,0,0,0];

%====H=HH'*HH===========
%[X,D]=eig(H);
% HH=(X*sqrt(D))';%结果一致
%=======================

QH1=[0.5,0;0,1];QH2=[2,2;-1,1];
b=[0,0,0,0,-1]';
c=[0,0,0,0,1]';
O4=zeros(4,4);O2=zeros(2,2);O3=zeros(3,3);
A11=[O4,HH(:,1);HH(:,1)',0];
A12=[O2,QH1(:,1);QH1(:,1)',-q1(1)];
A1=-[A11,zeros(5,3),zeros(5,3);
    zeros(3,5),A12,zeros(3,3);
    zeros(3,5),zeros(3,3),O3];

A21=[O4,HH(:,2);HH(:,2)',0];
A22=[O2,QH1(:,2);QH1(:,2)',-q1(2)];
A2=-[A21,zeros(5,3),zeros(5,3);
    zeros(3,5),A22,zeros(3,3);
    zeros(3,5),zeros(3,3),O3];

A31=[O4,HH(:,3);HH(:,3)',0];
A32=[O2,QH2(:,1);QH2(:,1)',-q2(1)];
A3=-[A31,zeros(5,3),zeros(5,3);
    zeros(3,5),O3,zeros(3,3);
    zeros(3,5),zeros(3,3),A32];

A41=[O4,HH(:,4);HH(:,4)',0];
A42=[O2,QH2(:,2);QH2(:,2)',-q2(2)];
A4=-[A41,zeros(5,3),zeros(5,3);
    zeros(3,5),O3,zeros(3,3);
    zeros(3,5),zeros(3,3),A42];

A51=[O4,zeros(4,1);zeros(4,1)',1];
A5=-[A51,zeros(5,6);zeros(6,5),zeros(6,6)];

C1=[eye(4),zeros(4,1);zeros(1,4),0];
C2=[eye(2),zeros(2,1);zeros(1,2),-r1];
C3=[eye(2),zeros(2,1);zeros(1,2),-r2];
C=[C1,zeros(5,3),zeros(5,3);
    zeros(3,5),C2,zeros(3,3);
    zeros(3,5),zeros(3,3),C3];

%====H=HH'*HH===========
%[X,D]=eig(H)
% HH=(X*sqrt(D))'
%=======================
%[X,D]=eig(Q2)
% QQ2=(X*sqrt(D))'
%[X,D]=eig(Q1)
% QQ1=(X*sqrt(D))'
%=======================
format long
X02=[1,0,-0.5;0,1,0;-0.5,0,1];
X03=[180,0,-12;0,60,-2;-12,-2,1];
X0=[eye(5),zeros(5,3),zeros(5,3);
    zeros(3,5),X02,zeros(3,3);
    zeros(3,5),zeros(3,3),X03];
y0=[1;0;2;4;20];

[n,m]=size(A1);p=length(y0);
AM=zeros(n,m,p);
AM(:,:,1)=A1;
AM(:,:,2)=A2;
AM(:,:,3)=A3;
AM(:,:,4)=A4;
AM(:,:,5)=A5;

S0=C;
for i=1:5
    S0=S0-y0(i)*AM(:,:,i);
end

gama=0.9;epsl=0.00000001;
%S0=C-y0(1)*A1-y0(2)*A2-y0(3)*A3-y0(4)*A4-y0(5)*A5
%mdd=min(eig(S0))
end




% ==================================================
function [X,y,S,k] = sdp_pc(X0,y0,S0,Ag,b,C,gam,epsi)
%------------------------------------
global y_rec
%--------------------------------------
format long
disp(' ')
disp('Program sdp_pc.m')
% Data preparation.
b = b(:);
p = length(b);
n = size(C)*[1 0]';
n2 = n*(n+1)/2;
X = X0;
y = y0(:);
%-----------------------------------------------
y_rec=zeros(100,length(y0));
k_rec=1;
y_rec(k_rec,:)=y0';
%-----------------------------------------------
S = S0;
I = eye(n);
gap = sum(sum(X.*S))/n;
k = 0;
A = zeros(p,n2);
for i = 1:p
   A(i,:) = (svec(Ag(:,(i-1)*n+1:i*n)))';
end
% SDP iterations. 
while gap > epsi
  % Generate predictor direction.
   Xi = inv(chol(X));
   Si = inv(chol(S));
   E = kron_s(S,I);
   F = kron_s(X,I);
   x = svec(X);
   rp = b-A*x;
   rd = svec(C-S-mat_s(A'*y));
   Ei = inv(E);
   M1 = A*Ei;
   M = M1*F*A';
   Mi = inv(M);
   rc = svec(-0.5*(X*S+S*X));  
   f1 = F*rd - rc;
   dy = Mi*(rp+A*Ei*f1);
   ad = A'*dy;
   dx = -Ei*(f1 - F*ad);
   ds = rd - ad;
   DX = mat_s(dx);
   DS = mat_s(ds);   
   lx = min(eig(Xi'*DX*Xi));
   ls = min(eig(Si'*DS*Si));
     if lx >= 0
      al = 1;
     else 
      al = min([1 -gam/lx]);
     end
     if ls >= 0
      bl = 1;
     else
      bl = min([1 -gam/ls]);
     end
   Xp = X + al*DX;
   yp = y + bl*dy;
   Sp = S + bl*DS;
   sig = (sum(sum(Xp.*Sp))/sum(sum(X.*S)))^3;
   gap = sum(sum(X.*S))/n;
   tau = sig*gap;
   % Generate corrector direction.
   rc = svec(tau*I-0.5*(X*S+S*X+DX*DS+DS*DX));
   f1 = F*rd - rc;
   dy = Mi*(rp+A*Ei*f1);
   ad = A'*dy;
   dx = -Ei*(f1 - F*ad);
   ds = rd - ad;
   DX = mat_s(dx);
   DS = mat_s(ds);
   lx = min(eig(Xi'*DX*Xi));
   ls = min(eig(Si'*DS*Si));
     if lx >= 0
      al = 1;
     else 
      al = min([1 -gam/lx]);
     end
     if ls >= 0
      bl = 1;
     else
      bl = min([1 -gam/ls]);
     end
   % Update {X, y, S}.
   X = X + al*DX;
   y = y + bl*dy;
   S = S + bl*DS;
%-----------------------------------------------
k_rec=k_rec+1;
y_rec(k_rec,:)=y';
%-----------------------------------------------
   % Evaluate duality gap.
   gap = sum(sum(X.*S))/n;
   k = k + 1; 
end
%-----------------------------------------------------
y_rec=y_rec(1:k_rec,:);
%----------------------------------------------------------
end

