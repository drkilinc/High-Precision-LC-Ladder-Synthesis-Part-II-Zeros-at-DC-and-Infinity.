function [num,errorn]=num_Z0(p,Z0,k)
% This function computes the numerator polynomial of an
% immitance function: Z(p)=Z0+sum{k(1)/[p-p(i)] 
% where we assume that Z0=A(1)/B(1)which is provided as input.
% 
% Input: 
%-------- poles p(i) of the immitance function Z(p)
%         as a MatLab row vector p
%-------- Residues k(i) of poles p(i) 
% Output:
%-------- num; MatLab Row-Vector 
% which includes coefficients of numerator polynomial of an immitance function. 
% num=Sum{k(j)*product[p-p(i)]} which skips the term when j=i 
%  
%----- Step 1: Determine total numer of poles n 
%
n=length(p);
nn=n-1;
%
%----- Step 2: Generation of numerator polynomials:
% numerator polynomial=sum of 
% sum of
% {Z0*[p-p(1)].[p-p(2)]...(p-p(n)]; n the degree-full product
% +k(1)*[p=p(2)].[p-p(3)]..[p-p(n)];degree of (n-1); the term with p(1)is skipped.                
% +k(2)*[p-p(1)].[p-p(3)]..[p-p(j-1)].[p-p(j+1)]..[p-p(n)];degree of(n-1)-the term with p(2)is skipped 
% +.............................................
% +k(j)*[p-p(1)].[p-p(2)]..[p-p(j-1)].[p-p(j+1)]..[p-p(n)];degree of (n-1)-the term with p(j)is skipped.       
% +.............................................
% +k(n)[p-p(1)].[p-p(2)]...[p-p(n-1)];degree of (n-1)-the term with p(n)is
% skipped.
%       
% Note that we generate the numerator polynomial within 4 steps.
% In Step 2a, product polynomial pra of k(1)is evaluated.
% In Step 2b, product polynomial prb of k(j)is evaluated by skipping the term when i=j.
% In Step 2c, product polynomial prc of k(n)is evaluated.
% In Step 2d, denominator of Z0 is generated.
%------------------------------------------------------------------------
% 
% Step 2a: Generate the polynomial for the residue k(1)
pra=[1];
for i=2:n

  simpA=[1 -p(i)];
% pra is a polynomial vector of degree n-1; total number of enrees are n.
  pra=conv(pra,simpA);% This is an (n-1)th degree polynomial.
end
na=length(pra);
% store first polynomial onto firs row of A i.e. A(1,:)
for r=1:na
  A(1,r)=pra(r);
end
% Step 2a: Compute the product for 2<j<(n-1)

for j=2:nn
    prb1=[1];
        for i=1:j-1
        simpB=[1 -p(i)];
        prb1=conv(prb1,simpB);
        end
    % Skip j th term
    prb2=[1];
        for i=(j+1):n
            simpB1=[1 -p(i)];
            prb2=conv(prb2,simpB1);
        end
        prb=conv(prb1,prb2);
        nb=length(prb);
%
% Store j polynomials on to j th row of A; i.e. A(j,:)
for r=1:nb          
A(j,r)=prb(r);
end
%       
    end
% Step 2c: Compute the product term for j=n
prc=[1];
for i=1:nn
    simpC=[1 -p(i)];
    prc=conv(prc,simpC);
end
nc=length(prc);
% store n the polynomial onto n the row of A(n,:)
for r=1:nc   
A(n,r)=prc(r);
end

%
%------------------------------------------------------------------------
for i=1:n
    for j=1:n
        C(i,j)=k(i)*A(i,j);
    end
end
%
%-------- Step 4: Generate the numerator as a MatLab row vector.
for i=1:n
D(i)=0;%Perform the sum operation to compute numerator polynomial
end
for j=1:n
    for r=1:n
    D(j)=D(j)+C(r,j);
    end;% Here is the numerator polynomial of length n.
end
[denom,errord]=denominator(p);
prd_n=Z0*denom; % this is n th degree polynomial vector with length n+1
a(1)=prd_n(1);
for i=2:(n+1)
    a(i)=D(i-1)+prd_n(i);
end
%    
num=real(a);
errorn=imag(a);