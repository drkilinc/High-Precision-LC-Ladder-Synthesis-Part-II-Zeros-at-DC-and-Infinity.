% Main Program Case_Study.m
clc;clear ALL;close ALL
% Synthesize Z(p)=a(p)/b(p):
a=[0,0.513976968750998,1.157898233543669,1.213142248104650,0.694401166177243,0.188516938299314,0]
b=[1,4.198433713455139,7.913422823138356,9.268979838320304,7.284067053656422,3.683494822490302,1]
[ka,qa,Highpass_Elements,Lowpass_Elements]=GeneralSynthesis_Yarman(a,b)
[CType,CVal,LC]=InceAyarSentez(a,b)
CVal=general_synthesis(a,b)
aa=a;bb=b;% Save original vectors a(p) and b(p)
%--------------------------------------------------------------------------
% Step 1a: Find out if F(p) is positive real.
pa=roots(a)
pb=roots(b)
% Step 1b: Generate even part of F(p)
[A,B]=even_part(a,b)
% Selecet eps_zero to ignore the small terms of A(p).
eps_zero=1e-8
% Step 1c: Refine the coefficiemnts of a(p) and b(p) via parametric
% approach
[a,b,ndc]=Check_immitance(a,b)
eps_a=norm(aa-a)
eps_b=norm(bb-b)
%--------------------------------------------------------------------------
% Step 2:Extraction of highpass elements
% Step 2.1: Extraction of first highpass element
[kr,R,ar]=Highpass_Remainder(a,b);
k(1)=kr;
% Step 2.1a: Define Fr(p)=R(p)/ar(p) as a new function F(p)=a(p)/b(p)
clear a; clear b;
R(6)=0;a=R;b=ar;
% Step 2.1.b: Refine the coefficients of new F(p)=a(p)/b(p)
[a,b]=Laddercorrection_onF(eps_zero,a,b)
eps_a1=norm(R-a/a(1))
eps_b1=norm(ar-b/a(1))
%--------------------------------------------------------------------------
% Step 2.2: Extraction of the second highpass element.
[kr,R,ar]=Highpass_Remainder(a,b)
k(2)=kr;
% Step 2.2.a:Define Fr(p)=R(p)/ar(p) as a new function F(p)=a(p)/b(p) 
clear a; clear b;
R(5)=0;a=R;b=ar;
% Step 2.2.b: Refine the coefficients of new F(p)=a(p)/b(p)
[a,b]=Laddercorrection_onF(eps_zero,a,b)
eps_a2=norm(R/ar(1)-a)
eps_b2=norm(ar/ar(1)-b)
% Step 2.3: Extraction of the third highpass element
[kr,R,ar]=Highpass_Remainder(a,b)
k(3)=kr
% Step 2.3a: Define Fr(p)=R(p)/ar(p) as a new function F(p)=a(p)/b(p)
clear a; clear b;
R(4)=0;a=R,b=ar
% Step 2.3.b: Refine the coefficients of new F(p)=a(p)/b(p)
[a,b]=Laddercorrection_onF(eps_zero,a,b)
eps_a3=norm(R-a/a(1))
eps_b3=norm(ar-b/a(1))
%--------------------------------------------------------------------------
% Step 2.4: Extraction of the fourth highpass element.
[kr,R,ar]=Highpass_Remainder(a,b)
k(4)=kr;
% Step 2.4.a:Define Fr(p)=R(p)/ar(p) as a new function F(p)=a(p)/b(p) 
clear a; clear b;
R(3)=0;a=R;b=ar;
% Step 2.4.b: Refine the coefficients of new F(p)=a(p)/b(p)
[a,b]=Laddercorrection_onF(eps_zero,a,b)
eps_a4=norm(R/ar(1)-a)
eps_b4=norm(ar/ar(1)-b)
%--------------------------------------------------------------------------
% Step 2.5: Extraction of the fifth highpass element
[kr,R,ar]=Highpass_Remainder(a,b)
k(5)=kr
% Step 2.5a: Define Fr(p)=R(p)/ar(p) as a new function F(p)=a(p)/b(p)
clear a; clear b;
a=R;b=ar;
% Step 2.5.b: Refine the coefficients of new F(p)=a(p)/b(p)
[a,b]=Laddercorrection_onF(eps_zero,a,b)
eps_a5=norm(R-a/a(1))
eps_b5=norm(ar-b/a(1))
%--------------------------------------------------------------------------
% Step 3: Extraction of Lowpass elements
if abs(a(1))<eps_zero;
    q=LowpassLadder_Yarman(a,b),
    'a(1)=0, F(p)=a(p)/b(p) is a minimum function'
end
if abs(b(1))<eps_zero;q=LowpassLadder_Yarman(b,a),
    'b(1)=0, F(p)=a(p)/b(p) is not a minimum function'
end








