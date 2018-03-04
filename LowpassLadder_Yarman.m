function q=LowpassLadder_Yarman(a,b)
% This function synthesizes the given minimum function F(p)=a(p)/b(p) as a
% lowpass LC ladder. Therefore, a(p) and b(p) must be properly generated.
% Step 1: Check the immitance verify that it is minimum function
eps_zero=1e-8;na=length(a);nb=length(b);
% Step 1: Generate Exact Lowpass Ladder Function F(p)=a(p)/b(p)
aa=a;bb=b;
if nb>2 [aa,bb,ndc]=Check_immitance(a,b);end
% Step 2:Remove poles at infinity 
%nplus1=length(b);n=nplus1-1;
nb=length(b);
if nb>2
n_infinity=nb-ndc;
%
for i=1:n_infinity-2
    [Q,a1,b1,ndc]=ExtractTrZero_infinity(eps_zero,aa,bb);
    q(i)=Q;
    clear a;clear b
    aa=a1;bb=b1;
end
% Now we have F1(p)=b1(p)/a1(p)
% Compute last residue q(n)
q(nb-1)=b1(1)/a1(2);
% Compute terminating constant q(n+1)
q(nb)=b1(2)/a1(2);
end
if nb==2
    q(1)=1/a(2);q(2)=b(2)/a(2);
end
end