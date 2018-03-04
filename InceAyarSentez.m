function [CType,CVal,LC]=InceAyarSentez(a,b)

%clear
%clc
%close all
%--------------------------------------------------------------------------
% PR Minimum Function Generation to test the programs
%W=[]; ndc=0;c=[1 2 3 4 5 6 7 8 9 10 11 12 13 14];
%a0=0.8;nc=length(c);
%[aa,bb]=Minimum_Function(ndc,W,a0,c);
%[aaa,bbb,ndc,f]=Check_immitance(aa,bb,1);
%a=aaa;
%b=bbb;
%--------------- Begining of Yarman Synthesis -----------------------------
[a,b,ndc]=Check_immitance(a,b);
na=length(a);
nb=length(b);

if abs(na-nb)<2
	n=nb;
else 
	% terminate synthesis function
	error('nb is different than na. Your Function F(p)=a(p)/b(p)is not proper for our synthesis pakage');
end

% Step 1; Extraction of transmission zeros at zero

%[a,b,ndc]=Check_immitance(aa,bb);
	rza{1}=roots(a);
	rzb{1}=roots(b);
	
	f=1;	% 'imp'
	CType=[];
	CVal=[];
	node_count=0;
	node_next=1;
	node_gnd=0;
	clear LC;
	NumFormatLength=3;
	
for i=1:ndc
	[k(i),a,b,f]=ExtractPoleAtDC_Yarman(a,b,f);
	if (ndc>i)
		a(end)=0;
	end
	
	% circuit marks for schematic and analysis
	node_count = node_next;
	if f==1
		CType	=[CType,2]; % serial C
		CName='C';
		node_in=node_count;
		node_out=node_count+1;
		node_next=node_count+1;
	else
		CType	=[CType,7];	% shunt L
		CName='L';
		node_in=node_count;
		node_out=node_gnd;
		node_next=node_count;
	end
	CVal	=[CVal,1/k(i)];
	% circuit list
	cj=length(CType);
	LC1=[CName, int2str(cj), ' ', int2str(node_in)];
	LC2=[' ', int2str(node_out), ' ', num2str(CVal(cj),NumFormatLength)]; 
    LC{cj} = [LC1,LC2];

	rza{i+1}=roots(a);
	rzb{i+1}=roots(b);
end

% Step 2: Extract transmission zeros at infinity
if (n-ndc)>0
	%q=LowpassLadder_Yarman(a,b,f);
	[a,b,ndc,f]=Check_immitance(a,b,f); % Given a(p) and b(p) find ndc and recompute a(p) and b(p) pricesiley
	na=length(a);
	nb=length(b);

	n_infinity=nb-ndc;% Total number of transmission zeros at infinity
	% Extract all the transmission zeros at infinity within a loop
	for i=1:n_infinity-2
		[Q,a,b,f]=ExtractTrZero_infinity(a,b,f);
		q(i)=Q;

		% circuit marks for schematic and analysis
		node_count = node_next;
		if f==1
			CType	=[CType,1]; % serial L
			CName='L';
			node_in=node_count;
			node_out=node_count+1;
			node_next=node_count+1;
		else
			CType	=[CType,8];	% shunt C
			CName='C';
			node_in=node_count;
			node_out=node_gnd;
			node_next=node_count;
		end
		CVal	=[CVal,q(i)];
		% circuit list
		cj=length(CType);
		LC1=[CName, int2str(cj), ' ', int2str(node_in)];
		LC2=[' ', int2str(node_out), ' ', num2str(CVal(cj),NumFormatLength)]; 
		LC{cj} = [LC1,LC2];
	end
	
	% first of last two component
	q(nb-1)=b(1)/a(2);
		node_count = node_next;
		if f==0
			CType	=[CType,1]; % serial L
			CName='L';
			node_in=node_count;
			node_out=node_count+1;
			node_next=node_count+1;
		else
			CType	=[CType,8];	% shunt C
			CName='C';
			node_in=node_count;
			node_out=node_gnd;
			node_next=node_count;
		end
		CVal	=[CVal,q(nb-1)];
		% circuit list
		cj=length(CType);
		LC1=[CName, int2str(cj), ' ', int2str(node_in)];
		LC2=[' ', int2str(node_out), ' ', num2str(CVal(cj),NumFormatLength)]; 
		LC{cj} = [LC1,LC2];
		
	% the last one: terminating resistance
	q(nb)=b(2)/a(2);
		node_count = node_next;
		CType	=[CType,9];	% shunt R
		CVal	=[CVal,q(nb)];
		CName='R';
		node_in=node_count;
		node_out=node_gnd;
		node_next=node_count;
		% circuit list
		cj=length(CType);
		LC1=[CName, int2str(cj), ' ', int2str(node_in)];
		LC2=[' ', int2str(node_out), ' ', num2str(CVal(cj),NumFormatLength)]; 
		LC{cj} = [LC1,LC2];

end
%k,q
Plot_circuit2(CType,CVal);
[Pay,Payda]=scam2ali(LC');
return

%-----------------------------------------------------------------------
function LC=ListCircuit(CType,CVal)
node=1;
ngnd=0;
for i=1:length(CType)
	switch (CType(i))
        case 1,    z='L'; nnode=node+1;	n1=node; n2=nnode;	% ed='L in Serial branch                 ';
        case 2,    z='C'; nnode=node+1;	n1=node; n2=nnode;	% ed='C in Serial branch                 ';
        case 3,    z='R'; nnode=node+1;	n1=node; n2=nnode;	% ed='R in Serial branch                 ';
        case 4,    z='L';		% ed='L of paralel R&L&C in Serial branch';
			if (i>1)&&((CType(i-1)==5)||(CType(i-1)==6))
				
			else
				nnode=node+1;	n1=node; n2=nnode;	
			end
        case 5,    z='C'; 			% ed='C of paralel R&L&C in Serial branch';
			if (i>1)&&((CType(i-1)==4)||(CType(i-1)==6))
				
			else
				nnode=node+1;	n1=node; n2=nnode;	
			end
        case 6,    z='R'; 			% ed='R of paralel R&L&C in Serial branch';
			if (i>1)&&((CType(i-1)==5)||(CType(i-1)==4))
				
			else
				nnode=node+1;	n1=node; n2=nnode;	
			end
        case 7,    z='L'; nnode=node+0;	n1=node; n2=ngnd;	% ed='L in shunt branch                  ';
        case 8,    z='C'; nnode=node+0;	n1=node; n2=ngnd;	% ed='C in shunt branch                  ';
        case 9,    z='R'; nnode=node+0;	n1=node; n2=ngnd;	% ed='C in shunt branch                  ';
        case 10,   z='L'; nnode=node+0; n1=node; n2=n1+1;		% ed='L of Serial R&L&C in shunt branch  ';
        case 11,   z='C'; nnode=node+2; n1=node+1; n2=ngnd;		% ed='C of Serial R&L&C in shunt branch  ';
        %case 12,   z='R'; nnode=node+1; n1=node; n2=ngnd;		% ed='R of Serial R&L&C in shunt branch  ';
        otherwise  z='X'; nnode=node+1;	n1=-i; n2=ngnd;	% ed='Non-defined element                ';
    end
    LC{i} = [z num2str(i) ' ' num2str(n1) ' ' num2str(n2) ' ' num2str(CVal(i))]; 
	node=nnode;
end


function [a1,b1,ndc,f]=Check_immitance(aa,bb,f)
% Given F(p)=a(p)/b(p)as a minimum function
% F(p) is free of finite jW zeros but it may have zeros at dc
% Find zeros of transmissions at dc
% Re-compute a(p) and b(p) to make F(p) minimum reactance

swap=0;
if poly_degree(aa) < poly_degree(bb)
	a=aa;	% impedance
	b=bb;
else
	a=bb;	% admittance
	b=aa;
	swap=1;
end

na=length(a);
nb=length(b); 
% Check the leading and the last coefficient of a(p)if they are zero
	if a(1)==0
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
%     sign=-1;
%     for i=na:-1:1
%         sign=-sign;
%         a_(i)=sign*a(i);
% 	end
	a_=polarity(a);
	
% 	sign=-1;
%     for i=nb:-1:1
%         sign=-sign;
%         b_(i)=sign*b(i);
% 	end
	b_=polarity(b);
	
	Num_Even=(conv(a,b_)+conv(a_,b))/2;
    n_Even=length(Num_Even);
    % change small value to zero
	% calculate the limit value of the zero
	LimitForZero=1e-3;%sum(abs(Num_Even))/length(Num_Even)/1000;
    if Num_Even(n_Even)==0
		for i=1:n_Even
			if abs(Num_Even(i))<LimitForZero
				Num_Even(i)=0;
			end
		end
    end
    % find degree of zero at dc
    for i=1:n_Even
        if abs(Num_Even(n_Even+1-i))>LimitForZero
            ndc=(i-1)/2;
            a0_sq=Num_Even(n_Even+1-i);
        end
    end
    if abs(Num_Even(n_Even))>LimitForZero
        ndc=0;
        a0_sq=Num_Even(n_Even);
        Comment_3='No dc transmission zero; ndc=0';
	end
	% find denom
    BB=conv(b,b_);
    B1=clear_oddpower(BB);
    nB=length(B1);
    %
	A1=a0_sq*R_Num(ndc,0);% A is specified in p-domain
	%A1=clear_oddpower(Num_Even);
	nA=length(A1);
	if (abs(nB-nA)>0)
		A1=fullvector(nB,A1);% work with equal length vectors
	end
	%
    [a1,b1]=RtoZ(A1,B1);

	% func type inverted
	if swap==1
		if (f==1)	%	'imp'
			f=0;	%	'adm'
		else
			f=1;	%	'imp'
		end
	end


function [kr,R,ar,f]=ExtractPoleAtDC_Yarman(a,b,f)
% In general, extract a pole at DC from a given immitance function F(p)
% Given a(p)=[a(1) a(2) ... a(n-1)   0]; a(1) could be zero or non-zero
% Given b(p)=[b(1) b(2) ... b(n-1) b(n)]
% F(p)=a(p)/b(p), F1(p)=b(p)/a(p)=k01/p+R(p)/a(p)
% k0=b(n)/a(n-1)
% R(p)=b(p)-k01*a(p)

n=length(b);

LimitForZero=sum(abs(a))/length(a)/1e-6;
if a(n)<LimitForZero
	kr=b(n)/a(n-1);
	
	% Set a(p)/p=[a(2)a(3)..a(n-1)]
	%			= a(2)p^(n-2)+a(3)p^(n-3)+..+a(n-2)p+a(n-1)
	
	% a(p)/p shift a(p) to right
	a1=a(2:end-1);
	a1=fullvector(n,a1);
	Reminder=abs(b-kr*a1);
	if kr<0
		Attention='Numerical Failure: Residue is negative'
	end
	R=Reminder(1:end-1);
	ar=a(1:end-1);
	
	% func type inverted
	if (f==1)	%	'imp'
		f=0;	%	'adm'
	else
		f=1;	%	'imp'
	end
end


% function q=LowpassLadder_Yarman(aa,bb,f)
% % This function synthesizes the given minimum function F(p)=a(p)/b(p) as a
% % lossless ladder. Therefore, a(p) and b(p) must be properly generated.
% %
% % Step 1: Check the immitance verify that it is minimum function and
% % find total number of transmission zeros at dc i.e find ndc.
% % zeros of transmissions at w=0
% 
% [aaa,baa,ndc]=Check_immitance(aa,bb,f); % Given a(p) and b(p) find ndc and recompute a(p) and b(p) pricesiley
% if f==1
% 	a=aaa;	% impedance
% 	b=baa;
% else
% 	a=baa;	% admittance
% 	b=aaa;
% end
% 
% na=length(a);
% nb=length(b);
% 
% % Step 2; Extraction of transmission zeros at infinity
% n_infinity=nb-ndc;% Total number of transmission zeros at infinity
% % Extract all the transmission zeros at infinity within a loop
% for i=1:n_infinity-2
% 	% Step 2: Extract transmission zeros at infinity
% 	[Q,a,b,ndc]=ExtractTrZero_infinity(a,b);
% 	q(i)=Q;
% 
% 	% func type inverted
% 	if (f==1)	%	'imp'
% 		f=0;	%	'adm'
% 	else
% 		f=1;	%	'imp'
% 	end
% end
% q(nb-1)=b(1)/a(2);
% q(nb)=b(2)/a(2);
% 


function [Q,a1,b1,f]=ExtractTrZero_infinity(aa,bb,f)
% Given F(p)=aa(p)/bb(p) as a minimum function
% Extract a transmission zero at infinity
% Given minimum function F(p)=aa(p)/bb(p)
if aa(1)<1e-8
    if bb(1)>0
        %flip over the function
        a=bb;
        b=aa;
        Q=a(1)/b(2);%quation
        c=conv(b,[Q 0]);
        nc=length(c);
        a=fullvector(nc,a);
        Reminder=a-conv(b,[Q 0]);
%         nR=length(Reminder);
%         for i=1:nR-2
%           R(i)=Reminder(i+2);
% 		end
		R=Reminder(3:end);
%         clear b1
%         nb=length(b);
%         for i=1:nb-1
%             b1(i)=b(i+1);
%         end
%         clear a1
		b1=b(2:end);
        a1=R;
            
		% func type inverted
		if (f==1)	%	'imp'
			f=0;	%	'adm'
		else
			f=1;	%	'imp'
		end
	end
end
if aa(1)>0
    Comment='Some thing is wrong'
end
[a1,b1,ndc,f]=Check_immitance(a1,b1,f);







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
if ndc>0;
	a(na)=0;
end;


%--------------------
function n=poly_degree(p)
% exclude leading zeros and
% finds polynoms degree
LimitForZero= sum(abs(p))/length(p)/1e6;
nz=0; % number of leading zero
for i=1:length(p)
	if abs(p(i))>LimitForZero
		break;
	else
		nz=i;
	end
end
n=length(p)-nz;

function B=polarity(A)
% MatLab Programs for Parameteric Approach: Program 1: sign
% This program changes the sign of the coefficients for a given polynomial
% Polynomial A is given as an even polynomial in w domain.
% By setting p=-w^2 we compute the coefficients of the new even polynomial
% in p domain.
% Polynomial A is an even polynomial in w domain
% P(w^2)=A(1)[w^(2n)]+A(2)[w^2(n-1)]+A(3)[w^2(n-2)+...+A(n)[w]+A(n+1)
% Inputs
%-------------- MatLab Vector A=[A(1) A(2) A(3)...A(n) A(n+1)] :
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
 
function AT=clear_oddpower(AA)
% This function clears the odd power terms in a given MatLab Polynomial AA
na=length(AA);
r=fix(na/2);
        for j=1:r
	        AT(j)=AA(na-2*j+2);
        end
        for i=1:r
            ATT(r-i+1)=AT(i);
        end
        for i=1:r
            AT(i+1)=ATT(i);
        end
        AT(1)=AA(1); 

function B=fullvector(n,A)
na=length(A);
for i=1:na
    B(n-i+1)=A(na-i+1);
end
for i=1:(n-na-1)
    B(i)=0;
end



%

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

function Apoly=R_allzero(ndc,nz,W)

% This function computes the value of the RA
% RA(-p2)=(-1)^(ndc)*(p*p)^(ndc)*{[(p^2+W(1)^2]^2}...{[p^2+W(n)^2]^2}
% Inputs
%        ndc=transmission zeros at DC
%        nz=Transmission zeros at W(i)
%        W(i)=Transmission zeros at W(i)
%        p=jw a frequency point at which RA is computed.
% Output
%        A: Even polynomial coefficients of the numerator polynomial
%        R(-p^2)
%
% Computation Steps
%
% Initialization
Apoly=[1];
if (nz>0)
	for i=1:nz
		Wi2=W(i)*W(i);
		Cpoly=[1 0 Wi2];
		Apoly=conv(Apoly,conv(Cpoly,Cpoly));
	end
end
if(ndc>0)
	D=zeros(1,2*ndc+1);
	D(1)=(-1)^ndc;
else
	D=[1];
end
Apoly=conv(Apoly,D);



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
a=num_Z0(p,Z0,k);	% a=num;
b=real(poly(p));	% b=denom;


% 
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
prd(n)=prc;

function [num,errorn]=num_Z0_x(p,Z0,k)
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
% for i=1:n
% 	% exclude i.th term of p polynomial
% 	p_low=[];
% 	p_up=[];
% 	if i>1
% 		p_low=p(1:i-1).'; 
% 	end
% 	if i<n
% 		p_up=p(i+1:n).';
% 	end
% 	psub=[p_low,p_up];
% 	% evaluate polynom again by excluding root p(i)
% 	prr=poly(psub);
% 	AA(i,:)=prr;
% 	CC(i,:)=k(i)*prr;
% end
%------------
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

% Step 2b: Compute the product for 2<j<(n-1)
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
%max(abs(C-CC)./abs(CC))

%-------- Step 4: Generate the numerator as a MatLab row vector.
D=sum(C,1);		% numerator - Z0
% D=[];
% for i=1:n
% 	D(i)=0;%Perform the sum operation to compute numerator polynomial
% end
% for j=1:n
% 	for r=1:n
% 		D(j)=D(j)+C(r,j);
% 	end;% Here is the numerator polynomial of length n.
% end
%[denom,errord]=denominator(p);
denom=real(poly(p));
% prd_n=Z0*denom; % this is n th degree polynomial vector with length n+1
% a(1)=prd_n(1);
% for i=2:(n+1)
% 	a(i)=D(i-1)+prd_n(i);
% end
%
a=Z0*denom+fullvector(length(denom),D);

num=real(a);
errorn=imag(a);

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
%
%----- Step 2: Generation of numerator polynomials:
% numerator polynomial=sum of
% sum of
% {Z0*[p-p(1)].[p-p(2)]...(p-p(n)]; n the degree-full product
% +k(1)*[p-p(2)].[p-p(3)]..[p-p(n)];degree of (n-1); the term with p(1)is skipped.
% +k(2)*[p-p(1)].[p-p(3)]..[p-p(j-1)].[p-p(j+1)]..[p-p(n)];degree of(n-1)-the term with p(2)is skipped
% +.............................................
% +k(j)*[p-p(1)].[p-p(2)]..[p-p(j-1)].[p-p(j+1)]..[p-p(n)];degree of (n-1)-the term with p(j)is skipped.
% +.............................................
% +k(n)[p-p(1)].[p-p(2)]...[p-p(n-1)];degree of (n-1)-the term with p(n)is skipped.
%
for i=1:n
	% exclude i.th term of p polynomial
	p_low=[];
	p_up=[];
	if i>1
		p_low=p(1:i-1).'; 
	end
	if i<n
		p_up=p(i+1:n).';
	end
	psub=[p_low,p_up];
	% evaluate polynom again by excluding root p(i)
	pr=poly(psub);	% pr is a polynomial vector of degree n-1;
	C(i,:)=k(i)*pr;	% 
end
%-------- Step 4: Generate the numerator as a MatLab row vector.
D=sum(C,1);		% numerator - Z0
%
denom=real(poly(p)); % this is n th degree polynomial vector with length n+1
% add Z0 into numerator
a=Z0*denom+fullvector(length(denom),D);
% result
num=real(a);
errorn=imag(a);







