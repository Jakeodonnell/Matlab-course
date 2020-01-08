%target.m
function hits = targetB(x,y,n,a,b);
hits = 0;
for i = 1:n
    %if(((x(i)).^2 + (y(i)).^2) < 1)
    if((((x(i).^2)./(a.^2)) + (y(i).^2)./(b.^2)) <= 1)
        plot(x(i),y(i),'o')
        xlim([-a a]);
        ylim([-b b]);
        pause(0.1)
        hold on
        hits = hits + 1;
    end
end
hits = ((a*b) * hits/n);
end