function [fitparms,fitdata,chisq] = lsqfit113(fun,xdata,ydata,init,fix,err,lb,ub,options)
% LSQFIT113 fits data using lsqnonlin
%
%	FIT(FUN,XDATA,YDATA,INIT,FIX,ERR,LB,UB,OPTIONS) performs a nonlinear
%	least squares fit of the function FUN to XDATA and YDATA, with initial
%	values for the variable coefficients INIT. 
%   The Boolean vector FIX allows the user to fix certain parameters, by
%   setting the corresponding element to a nonzero value. ERR is a matrix
%   containing the statistical uncertainty associated with the
%   measurements, in terms of a standard deviation; it may be a 1-by-1
%   matrix, indicating that the uncertainty is fixed, or it may be 1xc,
%	rx1, or rxc, where r and c are the number of rows and columns in the
%	dependent data set. 
%   LB and UB are the lower and upper bounds place on each fit parameter;
%   bounds on fixed parameters are ignored. OPTIONS is
%	the structure used in the least squares algorithm.

if nargin < 3 %nargin = 'number of arguments in'
   error('Not enough inputs')
end

% Check the sizes of the arrays xdata and ydata to make sure they are
% arranged in column vectors

[rx,cx] = size(xdata);
[ry,cy] = size(ydata);

if rx ~= ry
   error('XDATA and YDATA do not have the same number of rows')
end

if rx == 1 & cx > 1
   error('Data must be in column vector format')
end

if (nargin < 4 | isempty(init))
   init = initfit(fun,xdata,ydata);
else
   if size(init,2)~=cy
      if size(init,2)==1
         init = repmat(init,1,cy);
      else
         error('Initialization array wrong size')
      end
   end
end
if (nargin < 5 | isempty(fix))
   fix = zeros(size(init));
elseif size(fix,2)~=cy
   if size(fix,2)==1
      fix = repmat(fix,1,cy);
   else
      error('Fixed parameter array wrong size')
   end
end

if (nargin < 6 | isempty(err))
   err = ones(ry,cy);
end
[rerr cerr] = size(err);
if rerr == 1
   err = repmat(err,ry,1);
end
if cerr == 1
   err = repmat(err,1,cy);
end
if size(err) ~= [ry cy]
   error('Uncertainty matrix wrong size')
end

if (nargin < 7 | isempty(lb))
   lb = -Inf*ones(size(init));
end
[rlb clb] = size(lb);
if clb == 1
   lb = repmat(lb,1,cy);
end
if (nargin < 8 | isempty(ub))
   ub = Inf*ones(size(init));
end
[rub cub] = size(ub);
if cub == 1
   ub = repmat(ub,1,cy);
end

if (nargin < 9 | isempty(options))
   %   options = optimset('Display','off');
   options = optimset;
end

for j = 1:cy
   m = 1;
   n = 1;
   fixparms = [];
   for k = 1:size(fix,1)
      if fix(k,j)
         fixparms(m) = init(k,j);
         m = m + 1;
      else
         varparms(n) = init(k,j);
         lbvar(n) = lb(k,j);
         ubvar(n) = ub(k,j);
         n = n + 1;
      end
   end
   if (cx == cy)
      fitvarparms(:,j)=lsqnonlin('fitfun113', varparms,lbvar,ubvar,options,fun,xdata(:,j),ydata(:,j),fix(:,j),fixparms,err(:,j)).';
   elseif (cx == 1)
      fitvarparms(:,j)=lsqnonlin('fitfun113', varparms,lbvar,ubvar,options,fun,xdata,ydata(:,j),fix(:,j),fixparms,err(:,j)).';
      else
      error('XDATA not the correct size')
   end
   
   m = 1;
   fitparms(:,j) = init(:,j);
   for k = 1:size(fix,1)
      if ~fix(k,j)
         fitparms(k,j) = fitvarparms(m,j);
         m = m + 1;
      end
   end
   if nargout >= 2
      if (cx == cy)
         fitdata(:,j) = fitfun113(fitparms(:,j),fun,xdata(:,j));
      elseif (cx==1)
         fitdata(:,j) = fitfun113(fitparms(:,j),fun,xdata);
      else
         error('XDATA wrong size')
      end
      

      if nargout == 3 
			nonanidx = find(~isnan(fitdata(:,j)));
         chisq(1,j) = sum(((ydata(nonanidx,j) - fitdata(nonanidx,j))./err(nonanidx,j)).^2);

      end
   end
end

function init = initfit(fun,xdata,ydata)
% INITFIT will calculate initialization parameters for the function fun
%

ry = size(ydata,1);
cy = size(ydata,2);

switch fun
case 'gaussian'
   [init(1,:) pos] = max(ydata);
   for m = 1:cy
      [tmp hmpos1] = min(abs(ydata(1:pos(m),m) - init(1,m)/2));
      if pos(m) < ry
         [tmp hmpos2] = min(abs(ydata(pos(m)+1:end,m) - init(1,m)/2));
         init(2,m) = xdata(hmpos2) - xdata(hmpos1);
      else
         init(2,m) = .5*max(xdata);
      end
   end
   init(3,:) = xdata(pos,1).';
case 'lorentzian'
   [init(1,:) pos] = max(ydata);
   for m = 1:cy
      [tmp hmpos1] = min(abs(ydata(1:pos(m),m) - init(1,m)/2));
      if pos(m) < ry
         [tmp hmpos2] = min(abs(ydata(pos(m)+1:end,m) - init(1,m)/2));
         init(2,m) = xdata(hmpos2) - xdata(hmpos1);
      else
         init(2,m) = .5*max(xdata);
      end
   end
   init(3,:) = xdata(pos,1).';
case 'zeroline'
   init(1,:) = ydata(ry,:)/xdata(ry);
otherwise
   error([fun,' must be initialized before fitting'])
end
