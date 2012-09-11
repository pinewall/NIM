% function [x,flag,res,iter,resvec]= Adapt_BST(A, b, tol, maxit, u, v, x0); 
%  Adaptive Block-Stiefel iterative method.
%
% Input parameter:
%   A, b: Given problem and its right hand side.
%
%   tol, maxit: stopping criteria. Iteration stop when 
%               norm(b-A*x)/norm(b)<= tol or 
%               the number of iterations > maxit. 
%   NOTE: the stoping condition is not checked. The parameters here
%       is simply for the uniform appearance of iterative methods.
%
%   u, v: the pair of parameters that used by the algorithm. They 
%       should be selected such that all the eigenvalues are contained 
%       in the interval [v, u] to guarantee the convergence. Although v
%       may be set so that it greater than some smalled eigenvalues, 
%       a good approximation is prefered. 
%
%   x0: Given initial guess of solution.
%
% Output parameter:
%   x: solution vector.
%   NOTE: since there is no residual nore check, we assume that the best
%       solution we could get is the last solution calculated when there is
%       no projection performed. When there is a projection, the residual
%       norm may increase. In this case, we will take which ever smaller
%       between projected and non-projected residual to determine the 
%       solution vector.
%
%   flag: flag = 0, given toterance is satisfied, 
%         flag = 1, given tol is not achieved.
%
%   relres: relative resudual norm(b - A*x)/norm(b) at exit.
%
%   iter: number of iterations is conducted. 
%
%   resvec: vector that contains the norm(b-Ax) for each iteration step.
%       NOTE: this is only for instrumentation.  
%   
% function [x,flag,res,iter,resvec]= Adapt_BST(A, b, tol, maxit, u, v, x0); 
%
function [x,flag,relres,iter,resvec]= Adapt_BST(A, b, tol, maxit, u, v, x0)
%
% check the input parameters and set default values
%
if (nargin < 6)
    disp('not enough input parameter');
    return;
end
n = size(A,1);
if (nargin < 7) 
    x0 = zeros(n,1);
end
r = b - A*x0;
normb = norm(b);
flag = 1;
iter = 0;
resvec = [];
%
% Initialization
%
alpha = 2.0/(u - v);
beta = (u + v)/(u-v);
gama = beta/alpha;
%
% the start-up iteration
%
while iter < maxit
    x = x0 + (1/gama) * r;
    r = b - A*x;
    del_v = x - x0;
    w = 2/gama;
    relres = norm(r)/normb;
    if relres<tol
        flag = 0;
        return
    end
    relres0 = relres;
    expected_iter = acosh(relres0/tol)/acosh(beta);
    %
    % iteration loop without check residual norm
    %
    % for j=1:min([maxit-iter, max(10,floor(expected_iter*0.6))])
    for j=1:min([maxit-iter, max(5,floor(expected_iter*0.5))])
    	w = 1/(gama - w/(4*alpha*alpha));
    	del_v = w*r + (gama*w - 1)* del_v;
    	x = x + del_v;
    	r = b - A*x;
    	if nargout == 5    % this block is only for instrumentation. no norm is needed
        	res = norm(r);
        	resvec = [resvec; res];
        	relres = res/normb;
    	end;
    end
    relres = norm(r)/normb;
    if relres < tol      % converged to the required tolerance
   	flag = 0;
   	iter = iter + j;
   	return;
    else                 % adjust the parameter v and repeat the iteration
       if relres <= relres0     % do it only if relres is reduced
    	beta = cosh(acosh(relres0/relres)/j);
	    v = (beta-1)*u/(beta+1);
        alpha = 2.0/(u - v);
	gama = beta/alpha;
       end
	x0 = x;
   	iter = iter + j;
    end;
end;
end   
