function [weights_f,weights_d,h,optpol] = optextrap(lambda,nfs,nds,varargin)
% function [weights_f,weights_d,h,optpol] = optextrap(lam,nfs,nds,varargin)
% Optimize the time step size iteratively over a convex minimax subproblem.
% The time stepper generating polynomials are extrapolated to guarantee
% accuracy.  The free threads allow degrees of freedom to optimize the
% resulting stability polynomial over a given curve \lambda.  The routine
% bisects over the time step size \h, attempting to minimize the maximum
% value of |R(\h*\lambda)|-1.  If the resulting minimax value is at most
% zero, the time step size has a feasible solution, so we set our lower
% bracket to \h and continue.  If there is no feasible minimax solution in
% the space of our polynomial generating function basis, we set the upper 
% bracket to \h and continue.  Each minimax problem is convex and reliably 
% solved with cvx.  Once our bisection brackets are within a specifiable
% tolerance we exit the loop and return the latest feasible solution.
nvar         = length(nfs);
ndep         = length(nds);
maxthread    = max(max(nfs),max(nds));
polyorder    = maxthread+1;

% Parse the varargin parameters
parsed       = optextrap_params(polyorder,lambda,varargin);
tol_bisect   = parsed.tol_bisect;
tol_feasible = parsed.tol_feasible;
h_min        = parsed.h_min;
h_max        = parsed.h_max;
max_steps    = parsed.max_steps;
do_plot      = parsed.do_plot;
verbose      = parsed.verbose;

P            = gbspolys(maxthread);
Pdep         = P(nds,:);
Pfree        = P(nfs,:);
Vdep         = rot90(vander(1./nds.^2));
Vfree        = rot90(vandern(1./nfs.^2,ndep));

% lam should be a column vector
if size(lambda,1) == 1
    lambda = lambda.';
end
% length(lam) should be at least nvar+1
assert(length(lambda)>=nvar+1,'Underdetermined: spectrum should contain at least nvar+1 values.')

% default return values
weights = [];
optpol = [];

if verbose>0
    fprintf(' [ init] h_min: %e h_max: %e\n',h_min,h_max);
end
for iter=1:max_steps
    % Stop bisecting when relative tolerance is achieved
    if ((h_max-h_min)/h_min < tol_bisect) || (h_max < tol_bisect)
        break;
    end
    h = (h_max + h_min)/2.;

    % Run the convex solvers
    optimizer = @(solver)least_deviation(h,lambda,Pfree,Pdep,Vfree,Vdep,solver);
    solvers = {'sdpt3','sedumi'};
    for ss=1:length(solvers)
        [success,optx,optval] = optimizer(solvers{ss});
        if ~success || optval==Inf || isnan(optval)
            if verbose>1; disp([solvers{ss} ' failed!']); end
        else
            break;
        end
    end

    % Check if CVX has found a feasible solution and bisect accordingly
    if success && optval<tol_feasible
        h_min = h;
        weights = optx;
    else
        h_max = h;
    end
    
    if verbose>0
        fprintf(' [%5.0d] h_min: %e h_max: %e\n',iter,h_min,h_max);
    end
 end

if isempty(weights)
    warning('Unable to find a feasible solution');
    weights_f = [];
    weights_d = [];
    return;
end

% Return the largest known feasible value
h = h_min;

% Compute the optimal extrapolated stability polynomial in high precision
Vd = rot90(vander(1./vpa(nds).^2));
Vf = rot90(vandern(1./vpa(nfs).^2,ndep));
b = vpa([1;zeros(ndep-1,1)]);
x = vpa(weights);
c = Vd\(b-Vf*x);
optpol = x.'*vpa(Pfree)+c.'*vpa(Pdep);
optpol = fliplr(optpol);

c = double(c);
optpol = double(optpol);

weights_f = weights;
weights_d = c;

% Plot the optimized stability contour
if do_plot
    figure;
    stability_contour(optpol);
end

end

function [success,weights,optval] = ...
            least_deviation(h,lambda,Pfree,Pdep,Vfree,Vdep,solver)
% function [status,weights,v] = ...
%            least_deviation(h,lambda,Pfree,Pdep,Vfree,Vdep,solver)
%
% Solve the least deviation problem \mininize \max|R(h\lambda)|
%   for specified set of \lambda, with order constraints: R(z)\approx \exp(z)
    precision='best';

    ncoefs = size(Pfree,2);
    pval = zeros(length(lambda),ncoefs);
    for i=1:length(lambda)
        pval(i,:) = (h*lambda(i)).^(0:ncoefs-1);
    end

    nvar = size(Pfree,1);
    ndep = size(Pdep,1);
    b = [1;zeros(ndep-1,1)];

    cvx_begin
        cvx_quiet(true)
        cvx_precision(precision)
        cvx_solver(solver)
        
        variable x(nvar)           % CVX will optimize over x, the free thread weights

        c = Vdep\(b-Vfree*x);      % Compute the dependent thread weights
        pol = x.'*Pfree+c.'*Pdep;  % Compute the extrapolated stability polynomial
        R = abs(pval*pol.')-1;     % Compute the magnitude along the contour
        
        minimize max(R)            % Optimize the minimax value of |R(\h\lambda)-1|
    cvx_end

    success = strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved');
    weights = x;
    optval  = cvx_optval;
end

function results = optextrap_params(s,lam,optional_params)
% function results = optextrap_params(s,lam,optional_params)
%
% Set default optional and parameter values
    ip = inputParser;
    ip.FunctionName = 'optextrap_params';

    ip.addParameter('tol_bisect',1.e-3,@isnumeric);
    ip.addParameter('tol_feasible',1.e-9,@isnumeric);
    ip.addParameter('h_min',0,@isnumeric);
    ip.addParameter('h_max',2.01*s^2*max(abs(lam)),@isnumeric);
    ip.addParameter('max_steps',1000,@isnumeric);
    ip.addParameter('do_plot',false,@islogical);
    ip.addParameter('verbose',1,@isnumeric);

    ip.parse(optional_params{:});
    results = ip.Results;
end
