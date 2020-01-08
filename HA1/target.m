function hits = target(x,y,n);
hits = 0;
for i = 1:n
    if(((x(i)).^2 + (y(i)).^2) < 1)
        %if((((x(i).^2)./(a.^2)) + (y(i).^2)./(b.^2)) <= 1)
        plot(x(i),y(i),'o')
        hold on
        hits = hits + 1;
    end
end
hits = 2*2*hits
end