%% Home assignment
%% Exercise 2.a
x = linspace(-2,2);
y = exp(-x.^2);
figure(1)
plot(x,y,'r')
hold on
grid on
y1 = exp(-(x-1).^2);
plot(x,y1,'b')
title('plot of exp(-x.^2) and exp(-(x-1).^2)')
xlabel('x')
ylabel('y')
legend( 'exp(-x.^2)',' exp(-(x-1).^2)')
%% Exercise 2.b
%a = x.^2 + x + 1;
%b = ((x-1).* (x-2));
%f = a./b;
x1 = linspace(-2,1);
x2 =linspace(1,2);
x3 = linspace(2,5);
y1 = ((x1.^2) + x1 + 1)./((x1-1).*(x1-2));
y2 = ((x2.^2) + x2 + 1)./((x2-1).*(x2-2));
y3 = ((x3.^2) + x3 + 1)./((x3-1).*(x3-2));
figure(2)
plot(x1,y1,'b',x2,y2,'b',x3,y3,'b')
axis([0 5 -100 80])


%% Exercise 3
t = 0:0.06:6;
tdot = 0:1:6;
Adot = [205 130 85 65 42 25 15];
A = 202*exp(-0.42*t);
figure(3)
%discrete
plot(tdot,Adot,'o','markers',4)
grid on
hold on
xlabel('time')
ylabel('A(t) (half-life)')
%function
plot(t,A,'r')

%halftime
f = @(t) (202*exp(-0.42*t))-101;
y1 = f(t);
plot(t,y1,'k')
halftime_at_time = fzero(f, 101)
legend( 'discrete','function', 'halftime function')

%% Exercise 4
r = linspace(0,50);
P = @(r) r.*exp(-(r/3)).*(1 - ((2*r)/3) + ((2.*(r.^2))/27));
f = @(r) (r.*exp(-(r/3)).*(1 - ((2*r)/3) + ((2.*(r.^2))/27))).^2;
int_p = integral(P,0,inf);
A = sqrt(1/int_p)
y = P(r);
plot(r,y), grid on
r1 = fzero(P, 0)
r2 = fzero(P,1.8182)
r3 = fzero(P,6.969)

%% Exercise 5
x = linspace(-2,2);
y = linspace(-2,2);
[X,Y] = meshgrid(x,y);
Z = exp(-(X.^2 + Y.^2));
figure(5)
mesh(X,Y,Z)
figure(6)
contour(X,Y,Z)

%% Exercise 6
f = @(x) x.^3;
x = linspace(0,1);
y = f(x);
for a = 1:15
    subplot(3,5,a)
    plot(x,y,'g')
    hold on;
    h = 10^(-a);
    df = (f(x+h)-f(x))./h;
    df2 = (f(x+h)-f(x-h))./(2*h);
    plot(x,df,'r',x,df2,'b')
    title("h = 10^" + -a)
end
%% Exercise 1a PART 2
%throwarrows.m
n = input('enter amount of throws');
x = -1 + 2*rand(1,n);
y = -1 + 2*rand(1,n);
targetCalc = target(x,y,n)
cirkelnsArea = pi

%% %% Exercise 1b
%throwarrows.m
n = input('enter amount of throws');
a = input('enter width of target');
b = input('enter height of target');
centers = [0 0];
radii = a;
axis square
viscircles(centers,radii);
x = -a + 2*a*rand(1,n);
y = -b + 2*b*rand(1,n);
targetCalc = targetB(x,y,n,a,b)
cirkelnsArea = (a./2)*(b./2)*pi

%% Exercise 2
%a
load('CCD.MAT');
%b
figure(1)
imagesc(C,[3,7])
colormap('gray')
%c
for i = 2:99
    for j = 2:99
        if((C(i,j) == 0))
            C(i,j) = median(median(C(-1+i:i+1,-1+j:1+j)));
        end
    end
end
figure(2)
imagesc(C,[3,7])
colormap('gray')

%% Exercise 3
tol = input(' Give value of tol ')
T0 = (650/4)*ones(100,100); % initial distribution T0
T0(1,:) = 100;
T0(100,:) = 200; % set edge temp
T0(:,1) = 100;
T0(:,100) = 250; % set edge temp
T1 = T0; % dimension T1
diff = Inf; % make sure to enter the loop
while diff > tol
    diff = 0; % set max difference to 0
    for i = 2:99
        for j = 2:99
            if (i > 49 && j > 29) && (i < 70 && j < 50)
            else
                T1(i,j)=(T0(i+1,j)+T0(i-1,j)+T0(i,j+1)+T0(i,j-1))/4;
                if abs(T1(i,j) - T0(i,j)) > diff
                    diff = T1(i,j) - T0(i,j); % Update difference
                end
            end
        end
    end
    imagesc(T1) % plot temp. distribution
    colormap('hot'), colorbar % color sclae
    pause(0.1)
    T0 = T1; % uppdate temp. distribution
end



%% Exercise 1 PART 3
x = -1:0.2:1;
y = [0.04 0.06 0.1 0.2 0.5 1 0.5 0.2 0.1 0.06 0.04];
xip = linspace(-1,1); % x-values between -1 and 1
yip = interp1(x,y,xip,'spline');
plot(x,y,'o',xip,yip)
%% Exercise 3
load('aktivitet.mat');
a0 = [5000 ; 0.05 ; 20000 ; 0.1];
format long
[a,n] = gaussnewton(@fun,@dfun,a0,t,y,1e-5)

%% Exercise 4 i
%a
load('co2.data');
%b
x = linspace(1981,1990,234);
y = co2;
plot(x,co2,'+')
hold on
grid on
%c
%i
t = linspace(1981,1990,234);
y = co2;
a = polyfit(x.',co2,2);
y = a(1) + a(2)*t + a(3)*(t.^2);
xm = linspace(1981,1990,234);
ym = polyval(a,xm);
plot(xm,ym)

%% Exercise 4 ii
%a
load('co2.data');
%b
x = linspace(1981,1990,234);
y = co2;
plot(x,co2,'+')
hold on
grid on
%c
%i
t = linspace(1981,1990,234);
y = co2;
a = polyfit(x.',co2,2);
y = a(1) + a(2)*t + a(3)*(t.^2);
xm = linspace(1981,1990,234);
ym = polyval(a,xm);
plot(xm,ym)


%% Exercise 4 iii
load('co2.data');
t = linspace(1,234,234);
plot(t,co2,'+')
hold on
grid on
t = t.';
y = co2;
x = (linspace(1,234,234)')
k = (18*pi)./234;
A = [x.^0 x.^1 x.^2 sin(k.*t) cos(k.*t) sin(2*k.*t) cos(2*k.*t)];
a = A\y
xm = linspace(1,234,234);
ym = a(1)*xm.^0 + a(2)*xm.^1 + a(3)*xm.^2 + a(4)*sin(k*xm)
+ a(5)*cos(k*xm) + a(6)*sin(2*k*xm) + a(7)*cos(2*k*xm);
plot(xm,ym)
%% Exercise 4 d)
load('co2.data');
x = [linspace(1,234,234)]';
y = co2;
k = (18*pi)./234;
A = [x.^0 x.^1 x.^2 sin(k*x) cos(k*x) sin(2*k*x) cos(2*k*x)];
a = A\y
xm = linspace(1,1000,1000);
ym = a(1)*xm.^0 + a(2)*xm.^1 + a(3)*xm.^2 + a(4)*sin(k*xm) + a(5)*cos(k*xm) + a(6)*sin(2*k*xm) + a(7)*cos(2*k*xm);
spl = interp1(xm,ym,xm,'spline');
plot(xm, spl)
t = 817;
st = spl(817)
y1 = a(1)*t.^0 + a(2)*t.^1 + a(3)*t.^2 + a(4)*sin(k*t) + a(5)*cos(k*t)
+ a(6)*sin(2*k*t) + a(7)*cos(2*k*t)

%% Exercise 5a

x = [0 25 50 75 100 125 150 175 200 225 250 275 300 325]';
y = [0 22 38 53 67 74 87 94 101 115 122 126 126 126]';
plot(x,y,'o',x,y)
hold on
xlabel('Height(cm)')
ylabel('Bounce(cm)')
%b
A = [x.^1 x.^0];
a = A\y
xm = linspace(0,325);
ym = a(1)*xm.^1 + a(2)*xm.^0;
plot(xm,ym)

%% Exercise 5c

t = [0 25 50 75 100 125 150 175 200 225 250 275 300 325]';
y = [0 22 38 53 67 74 87 94 101 115 122 126 126 126]';
a0 = [1 ; 0.005];
[a,n] = gaussnewton(@funB,@dfunB,a0,t,y,1e-5)
yplot = a(1)*t ./(1+a(2)*t);
%f
maxheight = (a(1)/a(2))
plot(t,y,'o',t,yplot)
xlabel('Release')
ylabel('Bounce')
title('Table tennis ball release')

