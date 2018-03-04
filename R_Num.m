function A=R_Num(ndc,W)
% This function generates the numerator polynomial A(-p^2)of R(-p^2)in
% p-domain.
% A(-p^2)=A(1)p^2m+A(2)p^2(m-1)+...+A(m)p^2+A(m+1) 
% Where m=ndc+nz
% Inputs:
% RA(-p^2)=(-1)^(ndc)*(p*p)^(ndc)*{[(p^2+W(1)^2]^2}...{[p^2+W(nz)^2]^2}
%
%       ndc: number of zeros at DC
%       W: Finite transmission zeros of R(-p^2)
% Output:
%       A(p^2)=[A(1) A(2)...A(m) A(m+1)]
nz=length(W);%if threre is no finite transmission zero set nz=0
if (W==[0]) nz=0;end 
% Generate Numerator polynomial of RA(-p^2) in p-domain
Apoly=R_allzero(ndc,nz,W);
% Clear odd indexed terms from Apoly
A=clear_oddpower(Apoly);
