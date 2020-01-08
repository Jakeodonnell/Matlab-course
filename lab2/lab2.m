%%
%uppg 1a
f = @(x) tan(x) - x + 1;
x = linspace(0,(3*pi));
y = f(x);
figure(1);
plot(x,y), grid on
x = fzero(f, 4.5)
x = fzero(f, 7.5)

%%
%uppg 1b
f = @(x) sin(x) - 0.3*exp(x);
x = linspace(1, 2);
y = f(x);
figure(2);
plot(x,y), grid on
%solve gave us a root at 0.5419, therefore we look close to 0.5
x = fzero(f, 0.5)
x = fzero(f, 1)

%%
%uppg 1c
f = @(x) (0.1*(x.^3) - 5*(x.^2) - x + exp(-1 * x));
x = linspace(50, 60);
y = f(x);
r = solve(f)
figure(3);
plot(x,y), grid on
x = fzero(f, 0.3)
x = fzero(f, 50)

%%
%uppg 2
% f = @(d) 1 - 4*(sqrt(d.^2 + 9) - 3) .* d - (sqrt(d.^2 + 9))

f = @(d) 4*(sqrt(d.^2 + 9) - 3) .* d - (sqrt(d.^2 + 9))

d = linspace(-10, 10);
y = f(d);
ylim([0, 103]);
figure(4);
plot(d,y)
grid on
fzero(f, 2)

%%
%uppg 3
R = 0.0820578;
T = 500;
n = 2;
P = 100;
a = 3.640;
b = 0.04267;

f = @(V) ((n*R*T)./(V-n*b)) - a*((n./V)).^2 - P;

V = linspace(-10, 10);

y = f(V);
ylim([-80, 80]);
figure(5);
plot(V,y)
grid on
fzero(f, 1)

%%
%uppg 4

r = 1;
pp = 700;
pv = 1000;
g = 9.81;


x = linspace(-1, 3);

%V = pi*x.*(x - ((x.^2)/3));
ballWeight = ((4*pi*r^3)/3) * pp;
f = @(x) pv*g*(pi*x.*(x - ((x.^2)/3))) - ballWeight;
y = f(x);
figure(6);
plot(x,y)
grid on

fzero(f, -0.3)
fzero(f, 0.3)
%%
% uppg 5

k = 1.2;
e = 0.8;
T0 = 625;
%T1;
Tair = 300;
h = 20;
o = 5.67*10^-8;
dX = 0.05;

x = linspace(0, 600);
f = @(x) (e*o).*(x.^4 - Tair.^4) + h.*(x-Tair) - (k/dX.*(T0-x));
y = f(x);
figure(7);
plot(x,y)
grid on
T1 = fzero(f, 450)

%% 
%uppg 5.(2)
syms x
f = x.^2 / (x-1);

dfx = diff(f,x,1)
diff(dfx,x)

%%
%uppg 6
syms x;
f = log(x);
taylor(f,x,1,'order',5) % order = n-1;

%%
%uppg 6.(2)

syms x;
f = @(x) -exp(-x) .* sin(4*x);
df = diff(f,x);
g = @(x) sin(4*x).*exp(-x) - 4*cos(4*x).*exp(-x);
frev = @(x) -1 * (-exp(-x) .* sin(4*x));

x = linspace(0,3);
[minlocal, yValLocal] = fminbnd(f,0,3)
[maxlocal, yValLocalMax] = fminbnd(frev,2,3)
[maxGlobal, yValGlobalMax] = fminbnd(frev,0,1.5)

y = f(x);
yderr = g(x);

figure(8)
plot(x,y)
grid on
hold on
plot(x,yderr,'--')


%DOTS--------------------------
%BLUE = MIN GLOBAL, RED = MIN LOCAL

%find x value at 
getGlobal_X = fzero(g, 0)
%dot at local minimaplot(minlocal,yValLocal,'o', 'MarkerFaceColor', 'r')
%dot at global minima
plot(getGlobal_X,f(getGlobal_X),'o', 'MarkerFaceColor', 'b')

%BLACK = MAX LOCAL, GREEN = MAX GLOBAL
plot(maxlocal,-1 * yValLocalMax,'o', 'MarkerFaceColor', 'k')
plot(maxGlobal,-1 * yValGlobalMax,'o', 'MarkerFaceColor', 'g')
%------------------------------
syms x;
f = @(x) cos(2*x) + sin(x) + exp(-x.^2);
df = diff(f,x)
g = @(x) cos(x) - (2*sin(2*x)) - (2*x).*exp(-x.^2);
frev = @(x) -1 * (cos(2*x) + sin(x) + exp(-x.^2));

x = linspace(-2,2);
[minlocal, yValLocal] = fminbnd(f,0,2)
[maxlocal, yValLocalmax] = fminbnd(frev,0,2)

y = f(x);
yderr = g(x);


figure(9);
plot(x,y)
grid on
hold on
plot(x,yderr, '--')
%DOTS--------------------------
%BLUE = GLOBAL, RED = LOCAL

%find x value at 
getGlobal_X = fzero(g, -1.5)
%dot at local minima
plot(minlocal,yValLocal,'o', 'MarkerFaceColor', 'r')
%dot at global minima
plot(getGlobal_X,f(getGlobal_X),'o', 'MarkerFaceColor', 'b')

%BLACK = MAX LOCAL, GREEN = MAX GLOBAL
plot(maxlocal,-1 * yValLocalmax,'o', 'MarkerFaceColor', 'k')
%------------------------------

syms x;
f = @(x) (1./((x-3).^2 + 1)) - (x./((x-1).^2 + 0.1));
df = diff(f,x)
g = @(x)(x.*(2*x - 2))./((x - 1).^2 + 1/10).^2 - 1./((x - 1).^2 + 1/10) - (2*x - 6)./((x - 3).^2 + 1).^2;
frev = @(x) -1 * ((1./((x-3).^2 + 1)) - (x./((x-1).^2 + 0.1)));

x = linspace(0,5);
[minlocal, yValLocal] = fminbnd(f,0,5)
[maxlocal,yValLocalMax] = fminbnd(frev,0,5)

y = f(x);
yderr = g(x);


figure(10)
plot(x,y)
grid on
hold on
plot(x,yderr, '--')

%DOTS--------------------------
%BLUE = MIN GLOBAL, RED = MIN LOCAL
%dot at local minima
plot(minlocal,yValLocal,'o', 'MarkerFaceColor', 'r')
%BLACK = MAX LOCAL, GREEN = MAX GLOBAL
plot(maxlocal,-1 * yValLocalMax,'o', 'MarkerFaceColor', 'k')
%------------------------------

%%
%uppg 8 
syms x;
f = x.^2.*sin(x);
int(f,x)

%%
%uppg 9 
syms x;

%Analytical
f = x.*sin(x.^2);
I0 = int(f,x, 0, 10)

%Numerically
f = @(x) x.*sin(x.^2);
I1 = integral(f,0,10)

%%
%uppg 10
% x(t), y(t) = (t*cos(5t), sin(3t)), 0<= t <= 3
syms t;
t = 0:0.1:3;
x = (t.*cos(5*t));
y =  sin(3*t);
figure(11);
plot(x,y)

fx = @(t) (t.*cos(5*t));
fy = @(t) sin(3*t);

flenght = @(t) (sqrt((cos(5*t) - 5*t.*sin(5*t)).^2 + (3.*cos(3*t)).^2));
Length = integral(flenght,0,3)

%%
%uppg11

syms t;
t = 0:0.1:2*pi;
x = 8*cos(t);
y =  sin(t);
figure(12);
plot(x,y)

%-8*sin(x)
%cos(x)

flenght = @(t) sqrt((-8*sin(t)).^2 + (cos(t).^2));
Length = integral(flenght,0,3)
 
%%
%uppg12

syms x;
c1 = 30 / (60.^2);
f= @(x) c1.*x.^2;
df = @(x) sqrt(1 + ((2*c1*x).^2));



x = linspace(-70,70);

y = f(x);
figure(13);
plot(x,y)
grid on;
hold on;
plot(60, 30 ,'bo')
plot(-60, 30 ,'bo')

L = integral(df,-60,60)

%%
%13

L = 3;
B = 4;
syms x;

x = linspace(0,3);
f = @(x) 0.02*sin(4*pi*x);
y = f(x);
figure(15);
plot(x,y)

df = @(x) sqrt(1 + ((2*pi*cos(4*pi*x))/25).^2);

Lx = integral(df,0,3)
Ltot = Lx*4

%%
%14

c = @(t) cos((pi*t.^2)./2);
s = @(t) sin((pi*t.^2)./2);

C = integral(c,0,5)
S = integral(s,0,5)

%% DEL2
% 1a
x = linspace(-2,2);
y = linspace(-2,2);
[X,Y] = meshgrid(x,y);
Z = sin(X.^2 + Y.^2);
figure(16)
mesh(X,Y,Z)
figure(17)
contour(X,Y,Z)
xlabel('x');
ylabel('y');
zlabel('z');

%%
%1b

%x^2 + y^2 = r^2 ger intervall:
r = 0:0.2:2;
%2*pi is a whole period, pi gives half, 40 is resolution
theta = linspace(0, 2*pi, 40);
[R, THETA] = meshgrid(r, theta);

%I och med cirklens ekvation, beräknas X och Y med radianer. 
X = R.*cos(THETA);
Y = R.*sin(THETA);
Z = sin(X.^2 + Y.^2);
figure(18)
mesh(X,Y,Z)
figure(19)
contour(X,Y,Z)
xlabel('x');
ylabel('y');
zlabel('z');
%%
%2

syms x y;
f = sin(x.^2 + y.^2);
gradf = jacobian(f, [x,y])
pretty(gradf)

%%
%3a
x = linspace(-2, 2);
y = linspace(-2, 2);
[X, Y] = meshgrid(x,y)
Z = sin(X.^2 + Y.^2);
[Fx, Fy] = gradient(Z, 0.1 , 0.1);
%countor kan ses som derivatan vid denna punkt dvs lutning av funktionen.
figure(20);
contour(Z)
hold on 
quiver(Fx, Fy)

%%
%3b
r = -2:0.2:2;
%2*pi is a whole period, pi gives half, 40 is resolution
theta = linspace(0, 2*pi, 40);
[R, THETA] = meshgrid(r, theta);

%I och med cirklens ekvation, beräknas X och Y med radianer. 
X = R.*cos(THETA);
Y = R.*sin(THETA);
Z = sin(X.^2 + Y.^2);
[Fx, Fy] = gradient(Z, 0.5 , 0.5);
figure(21)
contour(Z)
grid on
hold on 
quiver(Fx, Fy)

%%
%4a
r = 0:0.1:7;
theta = linspace(0, 2*pi, 40);
[R, THETA] = meshgrid(r, theta);
X = R.*cos(THETA);
Y = R.*sin(THETA);

Z = sqrt(49.00001 - (X.^2) - (Y.^2));
mesh(X,Y,Z)
hold on

syms x y z
f = sqrt(49 - (x.^2) - (y.^2));
tangentPlane = taylor(f, [x y z], [6 2 3], 'order', 2)
x = linspace(-10,10);
y = linspace(-10,10);
[X,Y] = meshgrid(x,y);
Z = (49/3 - (2*Y)/3 - 2*X);
xlim([-10 10])
ylim([-10 10])
zlim([0 10])
xlabel('x');
ylabel('y');
zlabel('z');
figure(22);
mesh(X,Y,Z)

scatter3(6, 2, 3)

%%
%4b

syms x y z
f = sqrt(49 - (x.^2) - (y.^2));
tangentPlane = taylor(f, [x y z], [6 2 3], 'order', 3)

%%
%5
syms x y
f = x - log(x - y.^2);
gradf= jacobian(f, [x y]);
[xst yst] = solve(gradf(1), gradf(2), x, y)

hessian= jacobian(gradf, [x y]);
hessian_st = double(subs(hessian, [x y], [1 0]));
eigs(hessian_st)
%%
%6
syms x
x = linspace (-4, 4);
f = @(x) sin(x);
y = f(x);
figure(23);
plot(x,y)
hold on

g = @(x) ((x.^2) + 2);
y1 = g(x);
plot(x,y1)

nint = 1000;
xint = linspace(-.5,.5,nint)';
y1int = f(xint);
y2int = g(xint);

d = sqrt((xint - xint.').^2 + (y1int - y2int.').^2);
d = d + diag(repmat(inf,nint,1));
[mind,ind] = min(d(:));
[i1,i2] = ind2sub([nint,nint],ind);
plot([xint(i1),xint(i2)],[y1int(i1),y2int(i2)],'-bo')

%%
%7
syms x y
f = (x+y).*exp(x);
I = int(int(f,x,0,2),y,-1,2)

f = @(x,y) (x+y).*exp(x);

I = dblquad(f, 0,2,-1,2)

%%
%8


x = linspace(0, 3);
y = linspace(0, 5);
[X, Y] = meshgrid(x,y);
Z = 20 + (30 ./ ((X-2).^2 + (Y-1).^2 + 1));

%a
figure(24);
mesh(X,Y,Z)
xlabel('x');
ylabel('y');
zlabel('z');

T = @(x,y) 20 + (30 ./ ((x-2).^2 + (y-1).^2 + 1));
%b

I = dblquad(T,0,3,0,5)./(15)
%%
%9

z = @(x,y) 1./(1+x+y);

Volume = dblquad(z,0,2,0,1)

%%
%10
I = dblquad(@integrand,-2,0,-1,1)
function z = integrand(x,y)

n = length(x);
z = zeros(1,n);

for i = 1:n
    if(x(i) >= -2 & x(i) <= 0) & (x(i)+1)^2 + y^2 < 1
       z(i) = ((x(i).^3).*y + (x(i).^2).*(y.^2));
    end
end
end