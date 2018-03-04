% Main Program Example1_synthesis.m
% This program checks the given minimum function F(p)=a(p)/b(p)
% Program generates ndc and re-generates a(p) and b(p)
%
clear
clc
close all
a =[0    0.2641    0.6161    0.6473    0.4267    0.0472   0]
b =[1.0000    2.3326    6.0786   10.0776    7.6155    2.4686    0.2731]
[a1,b1,ndc]=Check_immitance(a,b)
eps_a=a-a1
eps_b=b-b1
