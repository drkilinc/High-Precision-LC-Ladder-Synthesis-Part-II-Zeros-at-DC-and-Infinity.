function [kinf,R,br]=removepole_atinfinity(a,b)
% a(1)>0,b(1)=0,b(2)>0
eps_zero=1e-8;
if a(1)>eps_zero
    if abs(b(1))<eps_zero
kinf=a(1)/b(2);
V=conv(b,[kinf 0]);
nv=length(V);
Remainder=fullvector(nv,a)-conv(b,[kinf 0]);
nR=length(Remainder);
for i=1:nR-2
    R(i)=Remainder(i+2);
end
    for i=1:nR-2
        br(i)=b(i+1);
    end
  
    end
else 'F(p)=a(p)/b(p) does not have a pole at infinity'
end
end
        
