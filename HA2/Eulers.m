function yy = euler(f,a,b,y0,h,n,t0)
t(1)
for i = 1:n
    if(i == 1)
        yy(i) = y0 + h*f(t0,y(0));
    else
        yy(i) = y0 + h*f(t0,y(i-1));
    end
end

end
