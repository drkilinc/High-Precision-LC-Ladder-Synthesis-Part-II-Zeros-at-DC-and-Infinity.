function [a,b]=Laddercorrection_onF(eps_zero,a,b)
% Given F(p)=a(p)/b(p)which is suppose to describe a ladder LC network
% This function corrects the coefficients of an Immitance function F(p)=a/b
%
% Case 1: F(p)=a(p)/b(p) is not minimum
if abs(a(1))>eps_zero
    if abs(b(1))<eps_zero
        [kinf,a_min,b_min]=removepole_atinfinity(a,b);
        [a,b,ndc]=kinfFmin_ToF(kinf,a_min,b_min);
    end
end
%
% Case 2: F(p)=a(p)/b(p) is minimum
if abs(b(1))>0
    if abs(a(1))<eps_zero
        [a,b,ndc]=Check_immitance(a,b);
    end
end