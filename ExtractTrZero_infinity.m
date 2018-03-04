function [Q,a1,b1,ndc]=ExtractTrZero_infinity(eps_zero,aa,bb)
% Given F(p)=aa(p)/bb(p) as a minimum function
% Extract a transmission zero at infinity
% Given minimum function F(p)=aa(p)/bb(p)
if abs(aa(1))<eps_zero
    if abs(bb(1))>eps_zero
        %flip over the function
        a=bb;
        b=aa;
        Q=a(1)/b(2);%quatient
        c=conv(b,[Q 0]);
        nc=length(c);
        a=fullvector(nc,a);
        Remainder=a-conv(b,[Q 0]);
        nR=length(Remainder);
        for i=1:nR-2
          R(i)=Remainder(i+2);
        end
        clear b1
        nb=length(b);
        for i=1:nb-1
            b1(i)=b(i+1);
        end
        clear a1
        a1=R;       
    end
end
if abs(aa(1))>eps_zero
    Comment='Some thing is wrong'
end
[a1,b1,ndc]=Check_immitance(a1,b1);
end