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
end