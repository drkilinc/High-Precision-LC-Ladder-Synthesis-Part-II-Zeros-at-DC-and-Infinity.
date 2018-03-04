% Main Program Case_Study2.m
clc;clear ALL;close ALL
% Synthesize Z(p)=a(p)/b(p):
a=[0,0.513976968750998,1.157898233543669,1.213142248104650,0.694401166177243,0.188516938299314,0]
b=[1,4.198433713455139,7.913422823138356,9.268979838320304,7.284067053656422,3.683494822490302,1]
[A,B]=even_part(a,b)
%[ka,qa,Highpass_Elements,Lowpass_Elements]=GeneralSynthesis_Yarman(a,b)
%[CType,CVal,LC]=InceAyarSentez(a,b)
%CVal=general_synthesis(a,b)
aa=a;bb=b;% Save original vectors a(p) and b(p)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Step 1:Extraction of the "Lowpass Element" as a shunt capacitor C1=kinf
% H(p)=1/F(p)=b(p)/a(p) is an admittance function which starts with a shunt
% capacitor C1
[kinf,R,br]=removepole_atinfinity(b,a)
q(1)=kinf;
R1=R;b1=br;
CVal=general_synthesis(b1,R1)
% Step 2: Extraction of "5-highpass elements" of Fr(p)=R(p)/ar(p)
clear a; clear b;
a=br;b=R;
for i=1:5
[kr,R,ar]=Highpass_Remainder(a,b)
k(i)=kr;
clear a; clear b;
a=R;b=ar;
end
k(6)=R/ar
Highpass_Components=1./k
% Termination conductance G=R/ar
G=R/ar
Resistance=1/G