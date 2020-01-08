%pingisboll
x = [0 25 50 75 100 125 150 175 200 225 250 275 300 325]';
y = [0 22 38 53 67 74 87 94 101 115 122 126 126 126]';
plot(x,y,'o',x,y), hold on
xlabel('Height of release')
ylabel('Bounce height')
% A = [x.^1 x.^0];
A = [x.^1 (1 + x.^1)]
a = A\y
% xm = 0:1:325;
% ym = a(1)*xm.^1 + a(2)*xm.^0;
% plot(xm,ym)

ym = (a(1)*x)/(1+a(2)*x);
plot(x,ym)