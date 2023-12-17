function err = fitfun113(varparms,fun,dataind,datadep,fixtest,fixparms,stddev)
%fitfun113 is used by FIT to perform nonlinear fits to data
%
%Syntax: fitfun113(varparms,fun,dataind,datadep,fixtest,fixparms)
%	returns the difference between datadep and the
%	function fun with variable parameters varparms, where
%	the independent variables are given in dataind.  If
%	datadep is left off, then the function returns an array
%	given by the function evaluated for varparms and dataind.
%  The variables fixtest and fixparms are available to pass
%  fixed values for selected parameters.  The length of
%  fixtest is equal the number of required inputs to fun,
%  and to the sum of the lengths of fixparms and varparms.
%  Elements of fixtest are nonzero for those parameters which
%  are to be fixed. The parameter lists varparms and fixparms
%  are ordered according to the order of the inputs to fun. 

if nargin < 3
   error('Not enough arguments')
elseif nargin < 4
   eval(['err = ',fun,'(varparms,dataind);'])
else
   if nargin == 4
      parms = varparms;
   elseif (nargin > 4) & (nargin < 8)
      if sum(~fixtest) ~= length(varparms);
         error('Inconsistent number of free parameters')
      end
      k = 1;
      m = 1;
      for j = 1:length(fixtest)
         if fixtest(j)
            parms(j) = fixparms(k);
            k = k + 1;
         else
            parms(j) = varparms(m);
            m = m + 1;
         end
      end
   else
      error('Too many inputs to fitfun113')
   end
   eval(['err = ',fun,'(parms,dataind,datadep,stddev);'])
end

function err=doubleexponential(parms,dataind,datadep,stddev)
% parms(1): ratio of the two exponentials' amplitudes.
% parms(2): first decay time (ps)
% parms(3): second decay time (ps)
% parms(4): overall amplitude

ideal = log(parms(4)*(exp(-dataind/parms(2))+parms(1)*exp(-dataind/parms(3)))/(1+parms(1)));
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end


function err = gaussian(parms,dataind,datadep,stddev)
%	GAUSSIAN fits a gaussian
%	parms(1) 	amplitude 
%	parms(2)	half-width at 1/e of maximum
%	parms(3)    center

ideal = parms(1)*exp(-((dataind-parms(3))/parms(2)).^2);
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end

function err=cosinefit(parms,dataind,datadep,stddev)

ideal = parms(1)+parms(2)*cos(parms(3)+dataind*parms(4));

if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end

function err = exponentialdrop(parms,dataind,datadep,stddev)
%	parms(1) 	time t=0
%	parms(2)	time constant of drop
%	parms(3)    amplitude of drop
%   parms(4)    long-time offset

t=dataind-parms(1);
j1=find(t<=0);
j2=find(t>0);
ideal(j1) = parms(3)+parms(4);
ideal(j2) = parms(3)*exp(-t(j2)./parms(2))+parms(4);

if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal')./stddev;
   err = diff;
end


function err = sin_decay(parms,dataind,datadep,stddev)
% parms(1): overall amplitude
% parms(2): exponential decay contstant
% parms(3): angular frequency
% parms(4): phase
% parms(5): y-offset

ideal = parms(1) .* exp(-dataind./parms(2)) .* sin(parms(3).*dataind + parms(4)) + parms(5);
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end

function err = sin_squared_decay(parms,dataind,datadep,stddev)
% parms(1): overall amplitude
% parms(2): exponential decay contstant
% parms(3): angular frequency
% parms(4): phase
% parms(5): y-offset

ideal = parms(1) .* exp(-dataind./parms(2)) .* (sin(parms(3).*dataind + parms(4))).^2 + parms(5);
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end


function err = sinc_squared(parms,dataind,datadep,stddev)
% parms(1): overall amplitude
% parms(2): angular frequency
% parms(3): phase
% parms(4): y-offset

V = sym(dataind);
ideal = double(parms(1) .* (sinc(parms(2) * V + parms(3))).^1 + parms(4));
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end


function err = exponential_rise(parms,dataind,datadep,stddev)
% parms(1): overall amplitude
% parms(2): x-offset
% parms(3): time_const
% parms(4): y-offset

ideal = parms(1) * (-1 + 2*exp(-(dataind-parms(2))/parms(3))) + parms(4);
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end

function err = sin_squared_rise(parms,dataind,datadep,stddev)
% parms(1): overall amplitude
% parms(2): exponential decay contstant
% parms(3): angular frequency
% parms(4): phase
% parms(5): y-offset
% parms(6): envelope x-offset
% parms(7): envelope scale

ideal =  (1 - exp(-(parms(7)*dataind+parms(6))./parms(2))) .* (parms(1) .* sin(parms(3).*dataind + parms(4))).^2 + parms(5);
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end


function err = exp_decay(parms,dataind,datadep,stddev)
% parms(1) = amp
% parms(2) = 

ideal =  parms(1) * exp(-2/parms(2)*(dataind)) + parms(4);
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end


function err = abs_sinc(parms,dataind,datadep,stddev)
% parms(1): overall amplitude
% parms(2): angular frequency
% parms(3): phase
% parms(4): y-offset

V = sym(dataind);
ideal = double(parms(1) .* abs(sinc(parms(2) * V + parms(3))) + parms(4));
if nargin < 3
   err = ideal;
else
   diff = (datadep-ideal)./stddev;
   err = diff;
end
