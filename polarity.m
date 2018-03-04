function B=polarity(A)
% MatLab Programs for Parameteric Approach: Program 1: sign
% This program changes the sign of the coefficients for a given polynomial
% Polynomial A is given as an even polynomial in w domain.
% By setting p=-w^2 we compute the coefficients of the new even polynomial
% in p domain.
% Polynomial A is an even polynomial in w domain
% P(w^2)=A(1)[w^(2n)]+A(2)[w^2(n-1)]+A(3)[w^2(n-2)+...+A(n)[w]+A(n+1)
% Inputs
%-------------- MatlLab Vector A=[A(1) A(2) A(3)...A(n) A(n+1)] :
%Coefficients of the even polynomial in w domain
%
% Outputs
%-------------- MatLab Vector B=[B(1) B(2) B(3) ...B(n) B(n+1)
%
sigma=-1;
NN=length(A);
n=NN-1;
for i=1:NN
    j=NN-i+1;
    B(j)=-sigma*A(j);
    sigma=-sigma;
end

