function [a,b]=Minimum_Function(ndc,W,a0,c)
% Generate analytic form of Fmin(p)=a(p)/b(p)
% Generate B(-p^2)
        C=[c 1];
        BB=Poly_Positive(C);% This positive polynomial is in w-domain
        B=polarity(BB);% Now, it is transferred to p-domain        
% Generate A(-p^2) of R(-p^2)=A(-p^2)/B(-p^2)
nB=length(B);
        A=(a0*a0)*R_Num(ndc,W);% A is specified in p-domain
nA=length(A);
if (abs(nB-nA)>0)
  A=fullvector(nB,A);% work with equal length vectors
end
        
% Generation of minimum immitance function using Bode or Parametric method
        [a,b]=RtoZ(A,B);% Here A and B are specified in p-domain
        na=length(a);
        if ndc>0;a(na)=0;end;
end