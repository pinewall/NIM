% 2 Sides Preconditioned Block-Stiefel Iterative Method with k iterations.
% The original problem to be solved:   Ax = b 
% To do the precondtion, the problem is transformed into:
%       inv(M')*A*inv(M)*y = inv(M')*b  where M'*M is the 
%       approximation of A, and M*x = y.
% NOTE1: matrix A is SPD to be success.
% NOTE2: There is no termination check in this method. The number of iteration
%        needed (k) is predetermined. 

% Input parameter:
%   A, b: Given problem and its right hand side.
%   u, v: the pair of parameter that used by the algorithm, 
%           They should be selected such that all the eigenvalues are
%           contained in the [v, u].
%   k: number of iteations required.
%   s: Dimension of projector. It is the estimation of the number of 
%      eigenvalues smaller than v.
%   M: preconditioner. M'*M is the approximation of A.
%   x0: Given initial guess of solution.
%
% Output parameter:
%   x: solution vector,
%   flag: identify if the projection is performed.
%         =0, no projection is performed
%         =1, projection is performed.
%   resvec: vector that contains the norm(b-Ax) for each iteration step.
%   R: the collection of the most recent s residual vectors.
%   
%   function [x, flag, R,resvec]= pbst_2sides(A, b, u, v, k, s, M, x0); 
%
function [x,flag, R, resvec]= pbst_2sides(A, b, u, v, k, s, M, x0); 
%
% check the input parameters and set default values
%
if (nargin < 6)
    disp('not enough input parameter');
    return;
end
n = size(A,1);
if (nargin < 8) 
    %x0 = rand(n,1);
    x0 = zeros(n,1);
end
if (nargin < 7)
    M = speye(n);
end
if nargout < 4
    resvec = zeros(1,1);
else
    resvec = zeros(k,1);
end
R = zeros(n,s);
s_count = 0;
R_start = 1;
flag = 0;


alpha = 2/(u - v);
beta = (u + v)/(u-v);
gama = beta/alpha;
y0 = M*x0;
r0 = M'\(b - A*x0);
y1 = y0 + (1/gama) * r0;
x = M\y1;
r1 = M'\(b - A*x);
nr1 = norm(r1);
del_v = y1 - y0;
w = 2/gama;
%
% start iteration
%
for j=1:k
    w = 1/(gama - w/(4*alpha*alpha));
    del_v = w*r1 + (gama*w - 1)* del_v;
    y1 = y1 + del_v;
    % applying precondition (next two lines)
    x = M\y1;
    r1 = M'\(b - A*x);
    %
    if nargout >= 4
        % collect residual info. 
        resvec(j) = norm(M'*r1);
    end
    % collect vector of Q
    if j >= k-s & j < k
        R(:,R_start) = r1;
        R_start = R_start+1;
        s_count = s_count + 1;
    end
end
nrk = norm(r1);
if (nrk > nr1/cheb(k, beta))
    % projection performed
    disp('projection');
    Q1 = orth(R);
    y0 = y1 + Q1*((Q1'*(M'\(A*(M\Q1))))\(Q1'*r1));
    x0 = M\y0;
    r1 = M'\(b - A*x0); 
    if (norm(r1)< nrk) 
        % accept the projection result
        x = x0;
        flag = 1;
        if nargout >= 4
            resvec(k) = norm(M'*r1);
        end
    end
end
    