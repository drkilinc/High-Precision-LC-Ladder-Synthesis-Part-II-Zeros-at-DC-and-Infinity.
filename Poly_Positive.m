function B=Poly_Positive(b)
% This function computes a positive polynomial in w domain
% B(w^2)=(1/2){b^2(w)+b^2(-w)}
% Input:
%       b=[b(1) b(2)...b(n) b(n+1)]
% Output:
%       B=[B(1)*w^2n+B(2)*w^2(n-1)...B(n)w^2 B(n+1)] 
 b2=conv(b,b);% Take the square of b(w)
 b2_=polarity(b2);
 C1=b2+b2_;
 C2=C1/2;
 B=clear_oddpower(C2);
 