function [x,max] = golden_search(f,a,b,tol)
%% Golden Section search algorithm to find a local maximum of a function on a
%compact real set - based on wikipedia webpage.

%INSPUTS:
%  -f:   real function where the domain is a subset of the reals
%  -a:   lower bound to search
%  -b:   upper bound to search
%  -tol: tolerance value

%OUTPUTS:
%  -max: the value that the function f attains at the maximum
%  -x:   the maximizer of f

gr = (sqrt(5) + 1)/2;  %golden rule

c = b - (b - a)/gr;
d = a + (b - a)/gr;

while abs(b-a)>tol
    %update bounds
    if f(c)<f(d)
       b = d;
   else
       a = c;
    end
    %recalculate c and d
    c = b - (b - a)/gr;
    d = a + (b - a)/gr;
end

x   = (a+b)/2;
max = f(x);
end

