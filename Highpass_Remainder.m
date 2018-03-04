function [kr,R,ar]=Highpass_Remainder(a,b)
% Computation of remainder for highpass element extraction
%  Here it is assumed that a(n)=0, a(1)>0,b(1)>0
% Given a minimum function F(p)=a(p)/b(p) such that a(p)=p*a^(p).
% Then we express H(p)=b(p)/a(p)=[1/F(p)]=kr/p+R(p)/ar(p) where R(p)=b(p)-kr*a^(p)
nb=length(b);
kr=b(nb)/a(nb-1);
for i=1:nb-1
    ar(i)=a(i);
end
ka=fullvector(nb,kr*ar);
for i=1:nb-1
    R(i)=b(i)-ka(i);
end
end

