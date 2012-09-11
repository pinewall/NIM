% Preconditioned Block-Stiefel iterative method.
% The original problem to be solved:   Ax = b 
% By doing the precondtion, the problem is transformed into:
%       inv(M1*M2)*A*x = inv(M1*M2)*b  where M1*M2 is the 
%       approximation of A.
% NOTE: matrix A and preconditioner are all SPD to be success.
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
%   s: The Number of most recent residual vectors that to be used for 
%       projection operation. s should be small compare to the allowed
%       number of iterations or the size of problem even when the number
%       eigenvalues smaller than v are unknown. 
%
%   k: number of iterations before projection operation is conduct. 
%   NOTE: iteration will terminate at min(k, maxit). This is the only
%       ternimation condition. When k>maxit, no prejection is
%       performed.When k<maxit, algorithm returns after the projection is
%       done.
%
%   M1, M2: preconditioner
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
%   flag: flag = 0, itertation success, 
%         flag = 1, iteration failed.
%       NOTE: since there is no condition check implemented, flag will not
%       be equal to 1.
%
%   res: norm of residual, norm(b - A*x).
%
%   iter: number of iterations if success. 
%       NOTE: Since there is no termination check, it will be equalto
%       min(k, maxit).
%
%   resvec: vector that contains the norm(b-Ax) for each iteration step.
%       NOTE: this is only for instrumentation.  
%   
% function [x,flag,res,iter,resvec]= PBST(A, b, tol, maxit, u, v, s, ...
%            M1, M2, x0); 
%
function [x,flag,res,iter,resvec]= PBST(A, b, tol, maxit, u, v, s, k, ...
            M1, M2, x0)
%
% check the input parameters and set default values
%
if (nargin < 6)
    disp('not enough input parameter');
    return;
end
n = size(A,1);
if (nargin < 11) 
    %x0 = rand(n,1);
    x0 = zeros(n,1);
end
if (nargin < 10)
    M2 = speye(n);
end
if (nargin < 9 )
    M1 = speye(n);
end
if (nargin < 8 )
    k = n/20;
end
if (nargin < 7 )
    s = 10;
end
%
% Initialization
%
normb = norm(b);
alpha = 2/(u - v);
beta = (u + v)/(u-v);
gama = beta/alpha;
iter = 0;
if nargout==5
    resvec = zeros(min( k, maxit ),1);
end
s_count = 0;
Q_start = 1;
Q = zeros(n, s);
%
% the start-up iteration
%
r0 = b - A*x0;
r1 = M2\(M1\r0); % preconditioned residual
x1 = x0 + (1/gama) * r1;
r0 = b - A*x1;
r1 = M2\(M1\r0); % preconditioned residual
del_v = x1 - x0;
w = 2/gama;
%
% iteration loop
%
for j=1:maxit
    w = 1/(gama - w/(4*alpha*alpha));
    del_v = w*r1 + (gama*w - 1)* del_v;
    x1 = x1 + del_v;
    %r1 = b1 - A1*x1;
    r0 = b - A*x1;
    r1 = M2\(M1\r0); % preconditioned residual
    if nargout == 5
        res = norm(r0);
        resvec(j) = res;
    end;
    % collect the most recent residual vectors before projection
    if j >= k-s & j < k
        Q(:,Q_start) = r1;
        Q_start = mod(Q_start,s)+1;
        s_count = s_count + 1;
    end
    % do projection
    if j == k
        % ready for projection
        res0 = norm(r0);
        disp('projection0');
        %Q = orth(Q);
        x0 = x1 + Q*((Q'*(M2\(M1\(A*Q))))\(Q'*r1)); % projection
        r0 = b - A*x0;
        res = norm(r0)
%         disp('projection1');
%         Q = orth(Q);
%         x0 = x1 + Q*((Q'*(M2\(M1\(A*Q))))\(Q'*r1)); % projection
%         r0 = b - A*x0;
%         res = norm(r0)
        if res < res0
            x = x0;
        else
            x = x1;
            res = res0;
        end
        if nargout ==5
            resvec(j) = res;
        end
        flag = 0;
        iter = j;
        return
    end        
end   
%end
flag = 1;
iter = j;
x = x1;
res = norm(b-A*x);