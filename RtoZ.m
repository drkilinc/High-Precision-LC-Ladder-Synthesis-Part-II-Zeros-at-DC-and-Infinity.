function [a,b]=RtoZ(A,B)
% This MatLab function generates a minimum function Z(p)=a(p)/b(p)
%      from its even part specified as R(p^2)=A(p^2)/B(p^2)
%      via Bode (or Parametric)approach
% Inputs: In p-domain, enter A(p) and B(p)
% A=[A(1) A(2) A(3)...A(n+1)]; for many practical cases we set A(1)=0.
% B=[B(1) B(2) B(3)...B(n+1)]
% Output:
%        Z(p)=a(p)/b(p) such that
%        a(p)=a(1)p^n+a(2)p^(n-1)+...+a(n)p+a(n+1)
%        b(p)=b(1)p^n+b(2)p^(n-1)+...+b(n)p+a(n+1)
% Generation of an immitance Function Z by means of Parametric Approach
% In parametric approach Z(p)=Z0+k(1)/[p-p(1)]+...+k(n)/[p-p(n)]
% R(p^2)=Even{Z(p)}=A(-p^2)/B(-p^2) where Z0=A(n+1)/B(n+1). 
%
% Given A(-p^2)>0
% Given B(-p^2)>0
% 
BP=polarity(B);%BP is in w-domain
AP=polarity(A);%AP is in w-domain
% Computational Steps
% Given A and B vectors. A(p) and B(p) vectors are in p-domain
% Compute poles p(1),p(2),...,p(n)and the residues k(i) at poles p(1),p(2),...,p(n)
[p,k]=residue_Z0(AP,BP);
%
% Compute numerator and denominator polynomials
Z0=abs(A(1)/B(1));
[num,errorn]=num_Z0(p,Z0,k);
[denom,errord]=denominator(p);
% 
a=num;
b=denom;