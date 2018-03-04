function [h]=Plot_Circuit2(varargin)
%
% Plot_Circuit1(CType, elm_val, h);
%
% inputs
% CType: component type vector 
% elm_val: Values vector
% h: figre handles
%
% Component types
%         case 1,    'L in Serial branch                 ';
%         case 2,    'C in Serial branch                 ';
%         case 3,    'R in Serial branch                 ';
%         case 4,    'L of paralel R&L&C in Serial branch';
%         case 5,    'C of paralel R&L&C in Serial branch';
%         case 6,    'R of paralel R&L&C in Serial branch';
%         case 7,    'L in shunt branch                  ';
%         case 8,    'C in shunt branch                  ';
%         case 9,    'R in shunt branch                  ';
%         case 10,   'L of Serial R&L&C in shunt branch  ';
%         case 11,   'C of Serial R&L&C in shunt branch  ';
%         case 12,   'R of Serial R&L&C in shunt branch  ';


	CType=varargin{1};
	CVal=varargin{2};
		
    %open figure
	if nargin==2
	    h=figure;
	else
		h=varargin{3};
	end
	
	if isempty(CType)
		return
	end
	if nargin>3
		tol=varargin{4};
	else
		tol=0.0001;
	end
	
    % starting location of overall figure, 
    xr0=0; yr0=20; 
	% input offset for circuit
	xd=5; yd=40;
	% draw
    [xx,cn]=draw_circuit2(h, CType, xr0+xd, yr0+yd, 0);

    % connection line
    X=[xr0 xr0+xd;  xr0 xr0+xd ]';
    Y=[yr0 yr0;     yr0+yd yr0+yd ]';
    
    ct=CType;
    tx=CVal;

    % set figure labels, borders, ...
    line(X,Y,'Color','k','LineWidth',2);
    title('Circuit Schematic');
    set(gca,'XGrid','off');
    set(gca,'YGrid','off');
    set(gca,'XTickLabel',{''});
    set(gca,'YTickLabel',{''});
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'box', 'on');

    % print component values
    pos=get(h,'Position');

    ts=length(tx); % how many label

    % label count for each row
    if ts>12
        cs=4;
    else
        cs=3;
    end
    rs=round(ts/cs+0.4999); % how many rows

    axis([-5 xx -2*(rs-1) yr0+yd+yd/2+2])

    etx=ig_selti(ct,tx,1,tol);
    rw= fix(xx/(cs))+1;       %  column width
    rh=5;           % row height
    row=-1;    col=0;
    for i=1:ts
        row=row+1;
        if row>=rs,
            col=col+1;
            row=0;
        end
        if tx(i)>0, vcl='k'; else vcl='r'; end
        text( 0+col*rw, 15-row*rh, etx(i),'HorizontalAlignment','left','FontSize',8,'color',vcl);
    end

%     % increase picture size if required
%     if (xx>300)||(cs>3),
%         pos(3)=pos(3)+250;
%     elseif xx>200
%         pos(3)=pos(3)+150;
%     end
    
    set(h,'Position',pos);


return

%-----------------------------------------------------------------------

function [x,ii]=draw_circuit2(hf,ctype,x0,y0,Num_offset)
% Draw circuit schematic
% Dr. Ali KILINÇ, 17-11-2003 
%
%
% inputs
% ctype: circuit components, given by the types defined in synthesis prog.
%	One resonant circuit components must be given in order of L-C-R.
% x0,y0: drawing start position
% Num_offset : offset for component index, used for display

wid=5; % component width
ls=17; % serial comp length
lp=40; % paralel comp length
lw=2; % line width
nc=length(ctype); % number of component
fs=8;%fix(1000/(1+nc*ls)+8); % font size

cl='k'; % line color
cL='r'; % inductance color
cC='b'; % cap. color
cR='k'; % resistor color

x=x0+wid;
y=y0;
y2=lp/3;

X=[x0 x; x0   x]'; 
Y=[y  y; y-lp y-lp]';
line(X,Y,'Color',cl,'LineWidth',lw);

ii=1;
%for i=1:nc,
while (ii<=nc),
    T='';
    switch ctype(ii)
        case 1,    %            z='#1-Serial L 
            A=draw_L( x,y,0,wid,ls,cL,lw);
            ts=['L' num2str(ii+Num_offset)]; tj='center';
            T='s';
        case 2,    %            z='#2-Serial C
            A=draw_C( x,y,0,wid,ls,cC,lw);
            ts=['C' num2str(ii+Num_offset)]; tj='center';
            T='s';
        case 3,    %            z='#3-Serial R
            A=draw_R( x,y,0,wid,ls,cR,lw);
            ts=['R' num2str(ii+Num_offset)]; tj='center';
            T='s';
            
		case 4,    %            z='#1-Serial L 
            A=draw_L( x,y+y2,0,wid,ls,cL,lw);
            ts=['L' num2str(ii+Num_offset)]; tj='center';
            T='s';
		case 5,    %            z='#2-Serial C
            A=draw_C( x,y-y2,0,wid,ls,cC,lw);
            ts=['C' num2str(ii+Num_offset)]; tj='center';
            T='s';
		case 6,    %            z='#3-Serial R
            A=draw_R( x,y,0,wid,ls,cR,lw);
            ts=['R' num2str(ii+Num_offset)]; tj='center';
            T='s';

		case 7,    %            z='#7-Shunt L 
            A=draw_L( x,y,01,wid,lp,cL,lw);
            ts=['L' num2str(ii+Num_offset)]; tj='left';
            T='p';
        case 8,    %            z='#8-Shunt C  
            A=draw_C( x,y,01,wid,lp,cC,lw);
            ts=['C' num2str(ii+Num_offset)]; tj='left';
            T='p';
        case 9,    %            z='#9-Shunt R
            A=draw_R( x,y,01,wid,lp,cR,lw);
            ts=['R' num2str(ii+Num_offset)]; tj='left';
            T='p';
		case 10,    %            z='#7-Shunt L 
            A=draw_L( x,y,01,wid,lp/3,cL,lw);
            ts=['L' num2str(ii+Num_offset)]; tj='left';
            T='p';
		case 11,    %            z='#8-Shunt C  
            A=draw_C( x,y-y2,01,wid,lp/3,cC,lw);
            ts=['C' num2str(ii+Num_offset)]; tj='left';
            T='p';
		case 12,    %            z='#9-Shunt R
            A=draw_R( x,y-2*y2,01,wid,lp/3,cR,lw);
            ts=['R' num2str(ii+Num_offset)]; tj='left';
            T='p';
        otherwise  %            z='non-defined element  
            A=[x,y];
            ts='';tj='center';
            T='p';
    end
    % write component name
    text(A(1),A(2),ts,'HorizontalAlignment',tj,'FontSize',fs);
    
    % for leaving free space after a component
    % check the type of next component 
    if ii<nc, 
        nx=ctype(ii+1); % next component exist
    else
        nx=-1;
    end
    if (nx>0)&(nx<7)
        Tnx='s';
    elseif (nx>6)&(nx<12)
        Tnx='p';
    else 
        Tnx='';
	end
	
    %elemanýn ana hat üzerindeki baðlantý çizgileri
    % draw wirings
    if T=='s',
		% seri kolda parlel rezonans devresi varsa
		% sonraki elemaný daha uzaða çiz
% 		if (ii<nc)&& ...
% 			((ctype(ii)==4)||(ctype(ii)==5)||(ctype(ii)==6)) && ...
% 			((ctype(ii+1)~=4)&&(ctype(ii+1)~=5)&&(ctype(ii+1)~=6))
% 			% uzun artýþ
			x1=x+ls+wid;
% 		else
% 	        x1=x+ls+wid/2;
% 		end
		% seri uç uzatma ve nötr hattý çizme
		if ii<nc
	        X=[x+ls x1; x    x1]'; 
		    Y=[y    y ; y-lp y-lp]';
		end
		% seri kolda parlel rezonans devresinin baðlantýlarý
		if (ctype(ii)==4), 
			X=[X, [x, x+ls; x, x+ls]];
			Y=[Y, [y, y;y+y2, y+y2]];
		elseif (ctype(ii)==5), 
			X=[X, [x, x+ls; x, x+ls]];
			Y=[Y, [y, y;y-y2, y-y2]];
		end
		
    else % T=='p'
		% sonraki elemaný çizme mesafesi
        if Tnx=='s',
            x1=x+wid;
        else
            x1=x+ls;
		end
		
		% sonraki eleman için seri uç uzatma ve nötr hattý çizme
		if (ii<nc),
			if (ctype(ii)~=ctype(ii+1))&& ...
				( ...
				  ( ...
					((ctype(ii)==4)||(ctype(ii)==5)) && ...
					((ctype(ii+1)==5)||(ctype(ii+1)==6)) ...
				  ) || ( ...
					((ctype(ii)==10)||(ctype(ii)==11)) && ...
					((ctype(ii+1)==11)||(ctype(ii+1)==12)) ...
				  )...
				)
				%seri rezonans devresi için
				X=[];
				Y=[];
			else
				X=[x x1; x   x1]';
				Y=[y  y; y-lp y-lp]';
			end
		end
		% paralel kolda seri rezonans devresinin baðlantýlarý
		if (ctype(ii)==10)&&(ii<nc)
			%LR devresi ise, Bu eleman L, ve C'nin yerini çiz
			if (ctype(ii+1)~=11), 
				X=[X, [x ; x  ]];
				Y=[Y, [y-y2 ;y-2*y2]];
			end
		elseif (ctype(ii)==11)
			% CR devresi
			if (ii==1)
				% sentezin ilk elemaný C ise L'nin yerini çiz
					X=[X, [x; x  ]];
					Y=[Y, [y;y-y2]];
			elseif (ii<nc)
				% ilk eleman rezonansýn L'si deðilse
				if (ctype(ii-1)~=10),
					X=[X, [x; x  ]];
					Y=[Y, [y;y-y2]];
				end
			end
			% LC devresi ise, R'nin yerini çiz
			if (ii>1)&&(ii<nc)
				if (ctype(ii-1)==10)&&(ctype(ii+1)~=12),
					X=[X, [x   ; x     ]];
					Y=[Y, [y-3*y2; y-2*y2]];
				end
			elseif (ii>1)&&(ii==nc)
				if (ctype(ii-1)==10),
					X=[X, [x   ; x     ]];
					Y=[Y, [y-3*y2; y-2*y2]];
				end
			end
		end
		
    end 
    
    % draw line
        line(X,Y,'Color',cl,'LineWidth',lw);

%     % terminate serial circiut
%     if ((ii==nc)&(ctype(end)=0)),
%         line(X,Y,'Color',cl,'LineWidth',lw);
%         line([x1 x1],[y y-lp],'Color',cl,'LineWidth',lw);
%     end
    
	% update cursor
	if (ii<nc)
		% paralel rezonans devresi için
		if (ctype(ii)~=ctype(ii+1))&& ...
			( ...
			  ( ...	
				((ctype(ii)==4)||(ctype(ii)==5)) && ...
				((ctype(ii+1)==5)||(ctype(ii+1)==6)) ...
			  ) || ( ...
				((ctype(ii)==10)||(ctype(ii)==11)) && ...
				((ctype(ii+1)==11)||(ctype(ii+1)==12)) ...
			  )...
			)
		%seri rezonans devresi için

			% artýþ yok
		else
			x2=x;
			x=x1;
		end
	else
			x2=x;
			x=x1;
	end
	%
	ii=ii+1;
end
    

return

%------------------------------------------------------------------------

function A=draw_C(x,y,dir,wid,len,c,lw)
% grafik ekrana kondansator seklinin
% yatay veya düsey eksende cizilmesi

wid=abs(wid); % genislik
len=abs(len); % uzunluk
dist=(wid/4); % yapraklar arasi aralik

% yatay veya düsey eksende cizim icin 
if dir==0,
    x0=x; % yatay
    y0=y;
    s=1;
else
    x0=y; % dusey
    y0=x;
    s=-1;
end

a=((len-dist)/2);
b=(wid/2);

x1=x0+s*a;
x2=x0+s*(len-a);
x3=x0+s*len;

X(:,1)=[x0 x1]';
Y(:,1)=[y0 y0]';

X(:,2)=[x2 x3]';
Y(:,2)=[y0 y0]';

X(:,3)=[x1 x1]';
Y(:,3)=[y0-b y0+b]';

X(:,4)=[x2 x2]';
Y(:,4)=[y0-b y0+b]';

% text loc
t1=((min(x0)+max(x3))/2);
t2=y0+wid/2;
to=wid/3;
% draw
if dir==0,
    line(X,Y,'Color',c,'LineWidth',lw);
    A=[t1 t2+to];
else
    line(Y,X,'Color',c,'LineWidth',lw);
    A=[t2+to t1];
end

%-----------------------------------------------------------------
function A=draw_L(x,y,dir,wid,len,c,lw)
% grafik ekrana bobin seklinin
% yatay veya düsey eksende cizilmesi

wid=abs(wid); % genislik
len=abs(len); % uzunluk
n=10; % daire uzerindeki adim sayýsý
tur=4; % bobin tur sayýsý

% yatay veya düsey eksende cizim icin 
if dir==0,
    x0=x; % yatay
    y0=y;
    s=1;
else
    x0=y;  % dusey
    y0=x;
    s=-1;
end

% tur capi
r=(wid/3);

% tur sayýsý verilen uzunluga sýðmýyorsa, azalt
while (tur*r*2>=(len-wid/2)), tur=tur-1; end

% kol uzunlugu
a=((len-tur*2*r)/2);

x1=x0+s*a;
X=[x0 x1]';
Y=[y0 y0]';

for i=n:-1:1,
    xa(i)=s*r*(1-cos(pi/n*i));
    ya(i)=y0+s*r*sin(pi/n*i);
end
for i=0:tur-1
    X=[X' x1+xa+2*s*r*i]';
    Y=[Y' ya]';
end

X=[X' x0+s*len]';
Y=[Y' y0]';

% text loc
t1=((min(X)+max(X))/2);
t2=y0+wid/2;
to=wid/3;
% draw
if dir==0,
    line(X,Y,'Color',c,'LineWidth',lw);
    A=[t1 t2+to];
else
    line(Y,X,'Color',c,'LineWidth',lw);
    A=[t2+to t1];
end

%------------------------------------------------

function A=draw_R(x,y,dir,wid,len,c,lw)
% grafik ekrana kondansator seklinin
% yatay veya düsey eksende cizilmesi

wid=abs(wid); % genislik
len=abs(len); % uzunluk

% yatay veya düsey eksende cizim icin 
if dir==0,
    x0=x; % yatay
    y0=y;
    s=1;
else
    x0=y; % dusey
    y0=x;
    s=-1;
end

% kutunun eni
a=((wid)/3);
% kutunun boyu
b=((len-1.5*wid)/2);

x1=x0+s*b;
x2=x0+s*(len-b);
x3=x0+s*len;

X(:,1)=[x0 x1]';
Y(:,1)=[y0 y0]';

X(:,2)=[x2 x3]';
Y(:,2)=[y0 y0]';

X(:,3)=[x1 x1]';
Y(:,3)=[y0-a y0+a]';

X(:,4)=[x2 x2]';
Y(:,4)=[y0-a y0+a]';

X(:,5)=[x1 x2]';
Y(:,5)=[y0-a y0-a]';

X(:,6)=[x1 x2]';
Y(:,6)=[y0+a y0+a]';

% text loc
t1=((min(x0)+max(x3))/2);
t2=y0+wid/2;
to=wid/3;
% draw
if dir==0,
    line(X,Y,'Color',c,'LineWidth',lw);
    A=[t1 t2+to];
else
    line(Y,X,'Color',c,'LineWidth',lw);
    A=[t2+to t1];
end



%-----------------------------------------------

function [out, out_dsc]=ig_selti(x,varargin)
% adding name and components value to write bottom of circuit 
% Sentez Sirasinda kullanilan eleman Tipleri ve konumu
% Ali Kilinç - 07*09*2001
n=length(x); z=[]; e=[]; ot=0;
if nargin >=2,
    e=varargin{1};
end
if nargin>2,
    ot=varargin{2};  
end
if nargin>3,
    tol=varargin{3};  
else
	tol=0.001;
end
% sayýlarý string/karakter dizisi olarak yazým formatý
ns=fix(log10(1/tol)+2);
if ns<3
	ns=3;
elseif ns>8
	ns=8;
end
frm=['%1.' num2str(ns-1) 'f'];

for i=1:n,
	%eleman deðerini K(kilo),M(Mega),m(mili),u(micro) ile yaz
	nsi=ns;
	if (abs(e(i))>=1e6) 
		kt=1e6;
		un='M';
	elseif (abs(e(i))>=1e3)
		kt=1000;
		un='K';
	elseif (abs(e(i))<1)&&(abs(e(i))>1e-3)
		kt=1e-3;
		un='m';
	elseif (abs(e(i))<1e-3)
		kt=1e-6;
		un='\mu';
	else
		kt=1;
		un='';
	end
	%ee=num2str(e(i)/kt,frm);
	ee=num2str(e(i)/kt,ns); 
	
    switch x(i)
        case 1,    z='L'; br='H';       % ed='L in Serial branch                 ';
        case 2,    z='C'; br='F';       % ed='C in Serial branch                 ';
        case 3,    z='R'; br='\Omega';  % ed='R in Serial branch                 ';
        case 4,    z='L'; br='H';       % ed='L of paralel R&L&C in Serial branch';
        case 5,    z='C'; br='F';       % ed='C of paralel R&L&C in Serial branch';
        case 6,    z='R'; br='\Omega';  % ed='R of paralel R&L&C in Serial branch';
        case 7,    z='L'; br='H';       % ed='L in shunt branch                  ';
        case 8,    z='C'; br='F';       % ed='C in shunt branch                  ';
        case 9,    z='R'; br='\Omega';  % ed='C in shunt branch                  ';
        case 10,   z='L'; br='H';       % ed='L of Serial R&L&C in shunt branch  ';
        case 11,   z='C'; br='F';       % ed='C of Serial R&L&C in shunt branch  ';
        case 12,   z='R'; br='\Omega';  % ed='R of Serial R&L&C in shunt branch  ';
        otherwise  z='X'; br='#';       % ed='Non-defined element                ';
    end
    out{i} = [z num2str(i) '= ' ee ' ' un '' br]; 
end

return
