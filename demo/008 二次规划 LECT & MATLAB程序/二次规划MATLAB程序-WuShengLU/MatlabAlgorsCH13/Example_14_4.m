function Example_14_4
close all;clear all;clc;
global x_rec
[F0,FF,c,L,epsl] = data_ex14_4;
%epsl=1e-7;
[xs,Fs,fs,k] = projective_sdp(c,FF,F0,epsl);
kNum=k
r=xs(1:2)
s=xs(3:4)
dist=norm(r-s,2)

figure(1)
ezplot('-0.25*tx.^2-ty.^2+0.5*tx+0.75')
hold on
ezplot('-(5*tx.^2+5*ty.^2+6*tx.*ty)/8+5.5*tx+6.5*ty-17.5')
axis([-1.5 4 -1.5 6])

xtext=1;ytext=4.0;
text(xtext,ytext,'S','FontSize',16)
xtext=0;ytext=0;
text(xtext,ytext,'R','FontSize',16)

plot(x_rec(:,1),x_rec(:,2),'r')
plot(x_rec(:,1),x_rec(:,2),'.r')
plot(x_rec(:,3),x_rec(:,4),'k')
plot(x_rec(:,3),x_rec(:,4),'.k')

plot([r(1),s(1)],[r(2),s(2)],'g')
%plot([y(1),y(3)],[y(2),y(4)],'g')

plot(x_rec(1,1),x_rec(1,2),'or')
plot(x_rec(1,3),x_rec(1,4),'ok')
%plot(y0(1),y0(2),'or')
%plot(y0(3),y0(4),'ok')

plot(r(1),r(2),'*r')
plot(s(1),s(2),'*r')
hold on
xtext=r(1)+0.08;ytext=r(2)-0.1;
text(xtext,ytext,'r^*','Color','m','FontSize',12)
xtext=s(1)+0.1;ytext=s(2)+0.1;
text(xtext,ytext,'s^*','Color','m','FontSize',12)

str=['Distance between R and S =',num2str(dist)]
title(str)
xlabel('x_1 / x_3');ylabel('x_2 / x_4');


end

function [F0,FF,c,L,epsl] = data_ex14_4
[AM,X0,y0,S0,C,b,gamma,epsl] = data_ex14_2;
p=length(y0);
[n,m]=size(AM(:,:,1));
F0=C;
F1=-AM(:,:,1);
F2=-AM(:,:,2);
F3=-AM(:,:,3);
F4=-AM(:,:,4);
F5=-AM(:,:,5);
FF=[F1 F2 F3 F4 F5];
c=[0,0,0,0,1]';
L=1;
epsl=5e-5;
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

% Program: projective_sdp.m
% Title: Projective algorithm for homogenized SDP problems.
% Description: Implements the projective algorithm (Algorithm 14.4).
%              It finds a vector x that minimizes c'*x 
%                  subject to: F(x) >= 0              
%              where F(x) = x(1)*F1 + x(2)*F2 + ... + x(n)*Fn + F0
%              is positive semidefinite.
% Theory: See Practical Optimization Sec. 14.6.3. 
% Input:      
%    c: p-dimensinal cost vector in Eq.(4.1).
%   FF: m x m*p matrix of the form FF = [F1 F2 ... Fp]
%   F0: m x m constant matrix constraint F(x) >= 0
% epsi: tolerance for the gap between f(xk) and its lower bound
% Output:    
%   xs: solution for the decision variables
%   fs: value of the objective function at x = xs 
%   Fs: constraint matrix F(x) at x = xs
%    k: number of iterations at convergence
% Example: 
% Solve the SDP problem
% minimize c'*x
% subject to: F0 + x1*F1 + x2*F2 + x3*F3 + x4*F4 >= 0
% where c = [1 1 1 1]', and matrices Fi are defined in
% Example 14.3.
% Solution:
% Execute the commands:
% c = [1 1 1 1]'
% [F0,FF] = data_ex14_3
% epsi = 1e-6
% [xs,Fs,fs,k] = projective_sdp(c,FF,F0,epsi)
% =================================================
function [xs,Fs,fs,k] = projective_sdp(c,FF,F0,epsi)
global x_rec
disp(' ')
disp('Program projective_sdp.m')
% Generate a strictly feasible initial point
[m,nn] = size(FF);
p = round(nn/m);
pz = p + 1;
mz = m + 1;
ct = [c(:); 0];
dt = [zeros(p,1); 1];
FFe = [FF F0];
FFt = zeros(mz,pz*mz);
FFt(mz,end) = 1;
for i = 1:pz,
    FFt(1:m,(i-1)*mz+1:i*mz-1) = FFe(:,(i-1)*m+1:i*m);
end
I = eye(mz);
Xki = I;
[x,Xkp] = proj(Xki,FFt,mz,pz);
vi = min(eig(Xkp));
while vi <= 0,
   Xw = Xki*Xkp - I;
   rou = max(abs(eig(Xw)));
   gk = 1/(1+rou);
   Xki = Xki - gk*Xw*Xki;
   [x,Xkp] = proj(Xki,FFt,mz,pz);
   vi = min(eig(Xkp));
end
k = 0;
Xk = Xkp;
xk = x;
%==========================================
k_rec=1;
x_rec=zeros(1000,p);
x_rec(k_rec,:)=xk(1:p)/xk(p+1);
%=============================================
fk = (ct'*xk)/(dt'*xk);
Xki = inv(Xk);
Xkp = Xk;
ve = min(eig(Xkp));
er = 1;
% Iteration begins
while er >= epsi,
   Wk = Xki*Xkp - I;
   if ve <= 0,
      Yk = Xkp - Xk;
   else
      rkd = trace(Wk*Wk) - 0.99^2;
      if rkd >= 0,
         Yk = Xkp - Xk;
      else
         [Fx,Fwi,xc,xd] = find_cd(ct,Xki,FFt,mz,pz);
         c2 = ct'*xc;
         d2 = dt'*xd;
         e2 = ct'*xd;
         ck = ct'*xk;
         dk = dt'*xk;
         wa = dk^2 + rkd*d2;
         wb = -(ck*dk+e2*rkd);
         wc = ck^2 + rkd*c2;
         fw = sqrt(wb^2 - wa*wc);
         fk = (-wb - fw)/wa;
         qw = zeros(pz,1);
         for i = 1:pz,
             qw(i) = trace(Fx(:,(i-1)*mz+1:i*mz));
         end
         pw = ct - fk*dt;
         pf = pw'*Fwi;
         qf = Fwi*qw;
         lams = -(pf*qw)/(pf*pw);
         xs = lams*pf' + qf;
         Xkpf = zeros(mz);
         for i = 1:pz,
             Xkpf = Xkpf + xs(i)*FFt(:,(i-1)*mz+1:i*mz);
         end
         Yk = Xkpf - Xk;
      end
   end
   rou = max(abs(eig(Xki*Yk)));
   gk = 1/(1+rou);
   Xki = Xki - gk*Xki*Yk*Xki;
   Xk = inv(Xki);
   [xk,Xkp] = proj(Xki,FFt,mz,pz);
   ve = min(eig(Xkp));
   vf = (ct'*xk)/(dt'*xk);
   k = k + 1;
   er = vf - fk;
%==========================================
   k_rec=k_rec+1;
   x_rec(k_rec,:)=xk(1:p)/xk(p+1);
%=============================================
end
% Generate solution point
xs = xk(1:p)/xk(p+1);
fs = c'*xs;
Fs = Xkp(1:m,1:m);
x_rec=x_rec(1:k_rec,:);
end

% Program: proj.m
% Description: This function is required by programs 
% projective_feasi.m (Algorithms 14.3) and projective_sdp  
% (Algorithm 14.4). It performs an orthogonal projection of the  
% input matrix Xk > 0 onto the subspace E = range(F)in the  
% metric <.,.> with respect to inv(Xk).
% Theory: See Practical Optimization Secs. 14.6.2 and 14.6.3.
% Input:
%    Xk: a positive matrix of size mz x mz
%   FFt: an (m+1) x (m+1)*(p+1) constraint matrix 
%        for argumented decision variables [x; tau]
% Output: 
%     x: a vectior of dimension nz such that 
%        Xkp = x(1)*F1 + x(2)*F2 + ... + x(n+1)*Fn+1   
%        minimizes the norm || Y - Xk || with respect to inv(Xk)
%   Xkp: Xkp = x(1)*F1 + x(2)*F2 + ... + x(n+1)*Fn+1
% =============================================================
function [x,Xkp] = proj(Xki,FFt,mz,pz)
Fx = Xki*FFt;
Fw = zeros(pz,pz);
for i = 1:pz,
   for j = 1:pz;
      if i > j,
         Fw(i,j) = Fw(j,i);
      else
         Fw(i,j) = trace(Fx(:,(i-1)*mz+1:i*mz)*Fx(:,(j-1)*mz+1:j*mz));
      end
   end
end
qw = zeros(pz,1);
for i = 1:pz,
   qw(i) = trace(Fx(:,(i-1)*mz+1:i*mz));
end
x = inv(Fw)*qw;
Xkp = zeros(mz,mz);
for i = 1:pz,
   Xkp = Xkp + x(i)*FFt(:,(i-1)*mz+1:i*mz);
end
end

% Program: find_cd.m
% Description: This function is required by projective_sdp.m
% (Algorithm 14.4). It computes matrices Ck and Dk defined
% in Eqs.(14.91a) and (14.91b).
% Input:      
%   Xki: inverse of matrix Xk
%    ct: (p+1)-dimensinal vector in Eq.(4.2)
%    FF: an mz x mz*pz matrix of the form FF = [F1 F2 ... Fn]
%    F0: constant term in constraint F(x) >= 0
% Output:
%    Ck and Dk: matrices defined in Eq.(4.19)
%    xc and xd: coefficient vectors to express Ck and Dk as linear
%               combination of matrices F~i's (see Problem 14.15).
% =============================================================
function [Fx,Fwi,xc,xd] = find_cd(ct,Xki,FFt,mz,pz)
Fx = Xki*FFt;
Fw = zeros(pz,pz);
for i = 1:pz,
   for j = 1:pz;
      if i > j,
         Fw(i,j) = Fw(j,i);
      else
         Fw(i,j) = trace(Fx(:,(i-1)*mz+1:i*mz)*Fx(:,(j-1)*mz+1:j*mz));
      end
   end
end
Fwi = inv(Fw);
xc = Fwi*ct;
xd = Fwi(:,end);
end
