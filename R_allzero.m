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
if (nz>=1)
        for i=1:nz
        Wi2=W(i)*W(i);
        Cpoly=[1 0 Wi2];
        Dpoly=conv(Cpoly,Cpoly);
        Apoly=conv(Apoly,Dpoly);
    end
end
if(ndc>0)
    D(2*ndc+1)=0.0;%in MatLab we shift all the terms by one.
    for i=1:2*ndc
        D(i)=0.0;
    end
    D(1)=(-1)^ndc;
end
        if ndc==0
            D=[1];
        end
        Apoly=conv(Apoly,D);