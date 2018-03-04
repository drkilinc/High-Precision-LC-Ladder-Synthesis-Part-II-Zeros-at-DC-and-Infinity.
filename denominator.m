function [denom,errord]=denominator(p)
% This function computes the denominator polynomial of an
% immitance function from the given poles.
% It should be noted that this form of the denominator is normalized with
% the leading coefficient b(1)=1: denom=(p-p(1))(p-p(2))....(p-p(n))
% Input: 
%-------- poles p(i) as a MatLab row vector p
%-------- Residues k(i) of poles p(i) 
% Output:
%-------- denom; MatLab Row-Vector 
% ---which includes coefficients of the denominator polynomial of an immitance function. 
% denom=product[p-p(i)] 
% 
% --------- Step 1: Determine n 
n=length(p);
% -------- Step 2: Form the product term.
 pr=[1]; %Define a simple polynomila pr=1.
for j=1:n 
     simple=[1 -p(j)];%this is the simple polynomial [p-p(j)]
% Generate multiplication of polynomials startting with 1.
      pr=conv(pr,simple);
end
denom=real(pr);
errord=imag(pr);