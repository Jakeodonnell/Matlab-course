function I = sumComp(f,a,b,c,d,n,m)
I = 0.0;
A = ((b-a)*(d-c)) / (m*n);

w = 2*ones(n+1,n+1);
w(2:n,2:n) = 2*w(2:n,2:n);
w(1,1) = 1;
w(1,n+1) = 1;
w(n+1,1) = 1;
w(n+1,n+1) = 1;

for i = 1:1:m+1
    for j = 1:1:n+1
        TC = w(i,j) * f(i,j);
        I = I + TC;
    end
end
I = ((A/4) * I)
