function [a1,b1,ndc]=Check_immitance(a,b)
% Given F(p)=a(p)/b(p)as a minimum function
% F(p) is free of finite jW zeros but it may have zeros at dc
% Find zeros of transmissions at dc
% Re-compute a(p) and b(p) to make F(p) minimum reactance
na=length(a);
nb=length(b);
% Check the leading and the last coefficient of a(p)if they are zero
eps_zero=1e-5;
if abs(a(1))<eps_zero
	if b(1)>0
		Comment_1='Degree of numerator is one less than denominator and F(p) is minimum function';
	end
	if a(na)==0
		Comment_2='There is at least one DC transmission zeros';
	end
end
%
if b(1)<0
	Comment_3='Given F(p) is not positive real';
end
% Computation of Numerator of Even Part
sign=-1;
for i=1:na
	sign=-sign;
	a_(na-i+1)=sign*a(na-i+1);
end
sign=-1;
for i=1:nb
	sign=-sign;
	b_(na-i+1)=sign*b(nb-i+1);
end
A1=conv(a,b_);A2=conv(a_,b);
Num_Even=(A1+A2)/2;
normA=norm(Num_Even);
n_Even=length(Num_Even);
%
if Num_Even(n_Even)==0
	for i=1:n_Even
		if abs(Num_Even(i))/normA<eps_zero
			Num_Even(i)=0;
		end
	end
end
%
a0_sq=1;
for i=1:n_Even
	if abs(Num_Even(n_Even+1-i))/normA>eps_zero
		ndc=(i-1)/2;
		a0_sq=Num_Even(n_Even+1-i);
	end
end
if Num_Even(n_Even)/normA>eps_zero
	ndc=0;
	a0_sq=Num_Even(n_Even);
	Comment_3='No dc transmission zero; ndc=0';
end
BB=conv(b,b_);
B1=clear_oddpower(BB);
nB=length(B1);
%
A1=a0_sq*R_Num(ndc,0);% A is specified in p-domain
nA=length(A1);
if (abs(nB-nA)>0)
	A1=fullvector(nB,A1);% work with equal length vectors
end
%
if na>1;[a1,b1]=RtoZ(A1,B1);end
if na==1;a1=a;b1=b;end
a1=abs(a1);b1=abs(b1);
end