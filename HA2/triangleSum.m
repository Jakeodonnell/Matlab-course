function I = triangleSum(f,T,x,y)
[m,n]  = size(T);
I = 0;
for i = 1:m
    xv = x(T(i,:));
    yv = y(T(i,:));
    cx = (1/3) * sum(xv);
    cy = (1/3) * sum(yv);
    matris = [xv';yv';1 1 1];
    
    A = (1/2) * abs(det(matris));
    I = I + f(cx,cy) * A;
end
end


