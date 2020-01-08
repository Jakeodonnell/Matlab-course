%   gaussnewton.m
function [a1,n] = gaussnewton(fun,dfun,a0,x,y,tol)
n = 0;
norm_da = 1e30;
while norm_da > tol & n < 100
    A = dfun(a0,x);
    f = fun(a0,x);
    da = A\(y - f);
    a1 = a0 + da;
    a0 = a1;
    norm_da = norm(da);
    n = n + 1;
end