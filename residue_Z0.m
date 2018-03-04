function [p,k]=residue_Z0(A,B)
% This function computes the residues k of a given real part R(w^2)
% R(w^2)=A(w^2)/B(w^2)
% It should be noted the R(w^2) is given in w domain as an even function of
% w rather than complex variable p=sigma+jw. 
% Inputs:
% ------- A(w^2); full coefficient numerator polynomial
% ------- B(w^2); full coefficient denominator polynomial
%Step 1: Find p domain versions of A and B
% That means we set w^2=-p^2
AA=polarity(A);
BB=polarity(B);
% Step 2: Find the LHP poles of R(p^2)
x=roots(BB);
p=-sqrt(x);
% Step 3: Compute the product terms
prd=product(p);
n=length(p);
% Step4: Generate residues
for j=1:n
    y=p(j)*p(j);
    Aval=polyval(AA,y);
    k(j)=(-1)^n*Aval/p(j)/B(1)/prd(j);
end 
