function prd=product(p)
% Input: 
%-------- poles p(i) as a MatLab row vector p
% Output:
%-------- product prd as a MatLab row vector
% prd(j)=p(j)=(p(j)-p1)(p-p2)...(p(j)-p(j-1)).(p(j)-p(j+1))...(p(j)-p(n) 
% which multiplies the terms [p(j)-p(i)] skipping the term when j=i.
% 
% This function computes the product of
% prd(j)=p(j)=(p(j)-p1)(p-p2)...(p(j)-p(j-1)).(p(j)-p(j+1))...(p(j)-p(n)
% The above product form is used to compute residues for parametric
% approach

n=length(p);
nn=n-1;
%
% Compute square of the poles
        for i=1:n
            ps(i)=p(i)*p(i);
        end
% Step 1: Compute the first product term for j=1
% pra=(p(1)^2-p(2)^2).(p(1)^2-p(3)^2...(p(1)^2-p(n)^2)
%
pra=1;
for i=2:n
    pra=pra*(ps(1)-ps(i));
end
prd(1)=pra;
%
% Step 2: Compute the product for 2<j<(n-1)skipping the j th the term.
% prb=(p(j)^2-p(1)^2)...((p(j)^2-p(j-1)^2).(p(j)^-p(j+1)^2)..(p(j)^2-p(n)^2)
%
for j=2:nn
    prb1=1;
        for i=1:j-1
        prb1=prb1*(ps(j)-ps(i));
        end
        prb2=1;
        for i=(j+1):n
     prb2=prb2*(ps(j)-ps(i));
        end
     prb=prb1*prb2;
     prd(j)=prb;
end
%
% Step 3: Compute the product term for j=n
% prc=(p(n)^2-p(1)^2).(p(n)^2-p(2)^2...(p(n)^2-p(n-1)^2)
%
prc=1;
for i=1:nn
    prc=prc*(ps(n)-ps(i));
end
if n>0;prd(n)=prc;end