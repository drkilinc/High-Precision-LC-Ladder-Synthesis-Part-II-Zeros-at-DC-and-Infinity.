function B=fullvector(n,A)
na=length(A);
for i=1:na
    B(n-i+1)=A(na-i+1);
end
for i=1:(n-na-1)
        B(i)=0;
end
