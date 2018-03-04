function [CVal,CType,Eleman,node,pay2,payda2]=Synthesis_LongDiv(pay,payda,z0,f0,repcount,spi,in_node,gr_node,tol)
% Verilen bir RLC empedans fonksiyonunun
% long division yöntemi ile terimlerine açar
% Ali Kilinç  14-11-2001, 1-8-2003, 18-8-2003
%
%function [CVal,CType,Eleman,node]=
%   Synthesis_LongDiv(pay,payda,z0,f0,repcount,spi,in_node,gr_node)
%
% OUTPUTS:
% CVal : extracted component values, used for circuit drawing
% CType: extracted component types listed by SETipi at below, used for circuit drawing
% Eleman: Output file format in Pspice format, used for circuit analyse 
% node: last node number
%
% INPUTS:
% pay: nominator
% payda: denominator
% z0=normalisation impedance
% f0=normalisation freq.
% repcount: how many element do you wont to synthesise ( 0: for to synthesize all components)
% spi: permission for extraction at zero frequency, (ex: set to 1)
%       (spi=0 disable, spi>0 enable) 
%       if spi is enabled, in firstly prog. remove poles at zero.
% in_node: input point number for pspice circuit, (ex: set to 1)
% gr_node: ground point number for pspice circuit (ex: set to 0)

% Sentez Sirasinda kullanilan eleman Tipleri ve konumu
% Ali Kilinç - 07*09*2001- 15*04*2007

SETipi=[ ...
    '#1-Serial L :';
    '#2-Serial C :';
    '#3-Serial R :';
    '#7-Shunt L  :';
    '#8-Shunt C  :';
    '#9-Shunt R  :'];

% SETip=[ ...
%     '#1-Seri Kolda L          :';
%     '#2-Seri Kolda C          :';
%     '#3-Seri Kolda R          :';
%     '#4-Seri Kolda Paralel L  :';
%     '#5-Seri Kolda Paralel C  :';
%     '#6-Seri Kolda Paralel R  :';
%     '#7-Paralel Kolda L       :';
%     '#8-Paralel Kolda C       :';
%     '#9-Paralel Kolda R       :';
%     '#10-Paralel Kolda Seri L :';
%     '#11-Paralel Kolda Seri C :';
%     '#12-Paralel Kolda Seri R :'];
%----------------------------------------------------------------------------------------
% main program starts here

% input parameter check

if z0<=0,
    z0=1;
end
if f0<=0,
    f0=1;
end
if repcount<=0,
    repcount=2*max(length(pay),length(payda))-1;
end
if spi~= 0,
   spi=1;
end
if (gr_node<0);
    gr_node=0;
end
if (in_node<=gr_node)
    in_node=gr_node+1;
end
if (tol < 1e-6) 
    tol=1e-6;
elseif (tol > 0.1)
    tol=0.1;
end
    

%local variables
%-----------------
emp=1;  % function is an empedance
sp=0;   % sp=0  remove pole at infinitiy
        % sp=1  remove pole at zero

neg_en = 1; % negative component enable

% active node
node=in_node;        

% extracted component count
cc=0;

repcount1=repcount;

%--------------------------
% main program loop
while repcount1 >= 1,
    pay1=p_decdeg1(pay);
    payda1=p_decdeg1(payda);
    
    if neg_en==1
        if (min(size(find(imag(pay)~=0)))>0) || (min(size(find(pay<0)))>0) ...
                || (min(size(find(imag(payda)~=0)))>0) || (min(size(find(payda<0)))>0)
            if cc>0
                disp(['Up to now ' int2str(cc),' Component extracted but now']);
            end
            disp('Input function has negative or complex value.');
            keyb= input('Do you want to continue (Y/N):','s');
            
            if (keyb=='y') || (keyb=='Y')
                neg_en=0;
            else
                break
            end
        end
    end

    % if there is a pole at zero, remove it first.
    pz=testPoleAtZero(pay1,payda1);    
     if (pz)&(spi),   % check the permission
        % x=1/s        
        % do a zero - pole transformation
        [ pay1 payda1]=PoleZeroTr(pay1,payda1);
        sp=1;
        pay1=p_decdeg1(pay1);
        payda1=p_decdeg1(payda1);
    end
    
    % check the degrees if numerator or denominator has zero degree, 
    % all components has extracted, synthesis has finished,
    % terminate synthesis
    n1=length(pay1);
    n2=length(payda1);
    if (n1==0)|(n2==0), 
        break; 
    end
    
    if n1<n2,
        % convert imp to adm or reverse
        if emp==1, emp=0; else emp=1; end
        p0=pay1;
        pay1=payda1;
        payda1=p0;
        % check the degrees
        n1=length(pay1);
        n2=length(payda1);
    end
    
    % divide polynoms and find the element to remove
    [q r]=deconv(pay1,payda1);    
    elm=q(1);
    
    % error
    if length(q)>2, 
        disp('This function is not a driving point impedance')
        break
    end  

    % end of synthesis
    if elm==0, 
        break; 
    end
    
    if n1 > n2,
        % remove a L or C element
        [pay,payda]=f_dec(pay1,payda1, [q(1) ,0], [0,1], tol);
        cc=cc+1;
        if (sp==1),
            if (emp==1)
                CVal(cc)=1/elm;
                CType(cc)=2;  % #2-Serial C      
                % new component types
                Name{cc}=['C' int2str(cc) ];
                N1(cc,1)=node;
                N2(cc,1)=node+1;
                E_val(cc,1)=1/elm;
                node=node+1;
            else 
                CVal(cc)=1/elm;
                CType(cc)=7;  % #7-Shunt L      
                % new component types
                Name{cc}=['L' int2str(cc) ];
                N1(cc,1)=node;
                N2(cc,1)=gr_node;
                E_val(cc,1)=1/elm;
            end;
        else
            if (emp==1)
                CVal(cc)=elm;
                CType(cc)=1;  % #1-Serial L      
                % new component types
                Name{cc}=['L' int2str(cc) ];
                N1(cc,1)=node;
                N2(cc,1)=node+1;
                E_val(cc,1)=elm;
                node=node+1;
            else 
                CVal(cc)=elm;
                CType(cc)=8;  % #8-Shunt C      
                % new component types
                Name{cc}=['C' int2str(cc) ];
                N1(cc,1)=node;
                N2(cc,1)=gr_node;
                E_val(cc,1)=elm;
            end;
        end
        repcount1=repcount1-1;
        
    elseif n1==n2,
        % remove a resistance
        %   disp('A resistance has found')
        [pay,payda]=f_dec(pay1,payda1, [elm], [1], tol);
        
        if (emp==1),
            %if (abs(elm)>1e-6)&&(abs(elm)<1e6)  
                cc=cc+1;
                CVal(cc)=elm;
                % new component types
                Name{cc}=['R' int2str(cc) ];
                N1(cc,1)=node;
                E_val(cc,1)=elm;
                if n1==1,
                    % if it is last resistor, set it as shunt
                    CType(cc)=9;  % #9-Shunt R
                    % new component types
                    N2(cc,1)=gr_node;
                else
                    CType(cc)=3;  % #3-Serial R
                    % new component types
                    N2(cc,1)=node+1;
                    node=node+1;
                end
            %end
        else
            %if (abs(elm)<1e6) && (abs(elm)>1e-6),
                cc=cc+1;
                CVal(cc)=1/elm;
                CType(cc)=9;  % #9-Shunt R
                % new component types
                Name{cc}=['R' int2str(cc) ];
                N1(cc,1)=node;
                N2(cc,1)=gr_node;
                E_val(cc,1)=1/elm;
            %end
        end        
    end
    
    % reverse of zero - pole transform
    if (pz)&(spi),
        pay1=p_decdeg1(pay);
        payda1=p_decdeg1(payda);
        % x=1/s    
        % zero - pole transform
        [ pay payda]=PoleZeroTr(pay1,payda1);
        sp=0;
    end
   
end % all

%------------------------------------
% synthesis has completed

% reverse of impedance - admitance transform
if emp==1,
    pay2=pay;
    payda2=payda;
else
    pay2=payda;
    payda2=pay;
end

% ----------------------------
% de-normalisation of impedance and frequency
% L0=z0/(2*pi*f0);
% C0=1/(z0*2*pi*f0);
L0=z0/(f0);
C0=1/(z0*f0);

% old type
for i=1:cc
    switch CType(i)
        case {1,4,7,10} % L
            CVal(i)=CVal(i)*L0;
        case {2,5,8,11} % C
            CVal(i)=CVal(i)*C0;
        otherwise % R
            CVal(i)=CVal(i)*z0;
    end
end

% new type comp.
for i=1:cc
    switch Name{i}(1)
        case {'L'} % L
            E_val(i)=E_val(i)*L0;
        case {'C'} % C
            E_val(i)=E_val(i)*C0;
        otherwise % R
            E_val(i)=E_val(i)*z0;
    end
end

% write elements  into a list
for i=1:cc
    Elem{i}=[ Name{i} '  ' int2str(N1(i)) '  ' int2str(N2(i)) '  ' num2str(E_val(i),'%8.3e') ];
end

Eleman=Elem';

return

%********************************************************************************
function po=p_decdeg1(pol,fark)
% this function removes leading zeros in polynomials
%
% polinomun en yuksek dereceli terimi sifir ise
% polinom derecesini bir azaltir.
%
% Ali Kilinc  %  25-Ekim-1997 % 07-Eylul-2001

if nargin < 2,
    fark = 1e-9; 
end

np=length(pol);
ma=max(abs(pol));

% if array is empty
if np==0,
    po=[];
    return
elseif np==1,
    if ma==0,
        po=[];
    else
        po=pol;
    end
    return
end

% if array is zero valued
if ma==0,
    po=[];
    return
end

%
if abs(pol(1))< ma*fark
    po=pol(2:np);
    po=p_decdeg1(po);  % recursivly decrease
else
    po=pol;
end

return

%***********************************************

function [pay,payda]=f_dec(pay1,payda1,pay2,payda2,tol)
%  iki rasyonel fonksiyonun farkini alir.
%  subtracts two rational function
%
%  f=f1-f2
%  f1=pay1/payda1
%  f2=pay2/payda2
%  f=pay/payda
%  
%  ALI KILINC - 29.3.1995 
%
% 11-08-2007 added
%  remove negative coefficients if less than relative tolerance
% 
payda=conv(payda1,payda2);
pay10=conv(pay1,payda2);
pay11=conv(pay2,payda1);

n=length(pay10)-length(pay11);
p=zeros(1,abs(n));

if n>0 ,
    pay11=[p pay11];
elseif 0>n ,
    pay10=[p pay10];
end

pay=pay10 - pay11;

be=pay./(abs(pay10) + abs(pay11) + eps);
%bei=find(abs(be)<eps);
bei=find(abs(be)<tol);
pay(bei)=0;

%*******************************************

function pz=testPoleAtZero(pay,payda)
% test if a function has a zero/pole at zero

pz=0; % set as no pole at zero

% degree
n1=length(pay);
n2=length(payda);

% break if there is no pole or zero, 
if (n1==0)|(n2==0), 
    return 
end

% max values
m1=max(pay);
m2=max(payda);

%
if abs(pay(n1))/abs(m1) < eps, pz=1; end
if abs(payda(n2))/abs(m2) < eps, pz=1; end

return

%****************************************************
function [p1, p2]=PoleZeroTr(pay,payda)
% if a function has a zero/pole at zero
% do a zero - pole transform as: x=1/s     

% degree
n1=length(pay);
n2=length(payda);

n=max(n1,n2);

p1 = fliplr([zeros(1,n-n1)    pay]);
p2 = fliplr([zeros(1,n-n2)  payda]);
return

% ******************************************************
function [pay,payda]=f_min(pay1,payda1);
%  [pay,payda]=f_min(pay1,payda1)
%
%  Verilen bir rasyonel fonksiyonun pay ve paydasýnda bulunan
%  ortak kökleri uzaklaþtýrarak sadeleþtirme yapar
%
%  ALÝ KILINÇ - 29.3.1995 
%  kökleri ayýrmada baðýl duyarlýlýk: epss=1e-6;

epss=1e-6;

pay1=p_decdeg1(pay1);	%
payda1=p_decdeg1(payda1);%

kat1=pay1(1);
kat2=payda1(1);
pay1(:)=pay1(:)./kat1;
payda1(:)=payda1(:)./kat2;

r1=roots(pay1); n1=length(r1);  r10=sort(r1);
r2=roots(payda1); n2=length(r2); r20=sort(r2);

i1=1;
while i1 <= n1,
    i2=1;
    while i2 <= n2,
        x1=abs( r10(i1)-r20(i2) );
        x2=abs(r10(i1)) + abs(r20(i2)) + eps;
        xx=x1/x2;
        if xx < epss,
            n1=n1-1;
            for i11=i1:n1, r10(i11)=r10(i11+1); end;    
            n2=n2-1;
            for i11=i2:n2, r20(i11)=r20(i11+1); end;    
            i1=0;
            i2=0;   
            %break
            return
        end
        i2=i2+1;
    end
    i1=i1+1;
end

for i1=1:n1,  r11(i1)=r10(i1); end;
for i1=1:n2,  r21(i1)=r20(i1); end;

pay=(kat1/kat2).*poly(r11);
payda=poly(r21);

return



