% function [y]=cheb(i, x);
% compute value of degree i Tchebyshev polynomial at x.
function [y] = cheb(i,x)
if abs(x)<= 1
    y = cos(i*acos(x));
else
    y = cosh(i*acosh(x));
end
