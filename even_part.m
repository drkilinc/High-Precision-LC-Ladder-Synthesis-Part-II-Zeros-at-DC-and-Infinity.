function [A,B]=even_part(a,b)
% Generate even part R(p)=A(p)/B(p) of a given immitance function F(p)=a(p)/b(p)
% Computation of Numerator of Even Part
na=length(a);
nb=length(b);
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
    Num_Even=(conv(a,b_)+conv(a_,b))/2;
      A=clear_oddpower(Num_Even);

    n_Even=length(Num_Even);
       BB=conv(b,b_);
    B=clear_oddpower(BB);
end