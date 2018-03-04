% Main Program Case_Study3.m
clc;clear all;close all
format long
% Generate minimum immitance function F(p)=a(p)/b(p) for Case Study 4:
% Inputs:
ndc=8;W=0;a0=1;
c=[1 0.2 -0.3 0.4 0.5 0.6 -0.7 0.8 0.9 1 .2 ...
    -.3 .4 0.5 0.8 0.6 -0.6 0.7 0.35 -0.45 ...
    0.75 0.65 1.0 -1 1];
[a,b]=Minimum_Function(ndc,W,a0,c);
% Synthesize F(p)=a(p)/b(p) as a lossless ladder in resistive termination
[k,q,Highpass_Elements,Lowpass_Elements]=GeneralSynthesis_Yarman(a,b)
%[CType,CVal,LC]=InceAyarSentez(a,b)