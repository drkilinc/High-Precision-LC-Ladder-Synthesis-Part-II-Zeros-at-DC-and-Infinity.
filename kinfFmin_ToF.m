function [a,b,ndc]=kinfFmin_ToF(kinf,a_min,b_min)
% Given kinf, a_min, b_min
% Generate a highprecision F(p) in ladder form
[a_min,b_min,ndc]=Check_immitance(a_min,b_min);
V=kinf*conv([1 0],b_min);
nv=length(V);
a=V+fullvector(nv,a_min);
b=[0 b_min];
end