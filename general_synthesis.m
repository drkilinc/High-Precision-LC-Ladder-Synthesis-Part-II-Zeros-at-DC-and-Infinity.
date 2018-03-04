function CVal=general_synthesis(pay,payda)
%pay=g+h;
%payda=g-h;
% Input parameters of function
	%pay = [ 0.858  0.763  1.75  0.653  0.461  0];
	% payda = [ 0.892  .793  2.21  1.03  1.08  0.121  0.0859];
	z0 = 1;	% set not to use normalization
	f0 = 1;	% set not to use normalization
	repcount = 0;	% synthesize all function
	spi = 1;	% include poles at zero to synthesis
	in_node = 1;	% define circuit input node
	gr_node = 0;	% define circuit ground node
	tol = 0.01;	% relative tolerance;
		% pay-payda coefficients have 3 digits

	% call the function
	[CVal,CType,Eleman,node,pay2,payda2]=Synthesis_LongDiv(pay,payda,z0,f0,repcount,spi,in_node,gr_node,tol)
Plot_circuit2(CType,CVal)