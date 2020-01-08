%% uppg1a
load freefall.mat

a = 1;
b = 200;
n = 201;
h = (b-a)/(n-1);
t = a:h:b;

dy(1) = (y(2)- y(1))/h;
dy(201) = (y(201) - y(200))/h;

for i=2:n-1
    dy(i) = (y(i+1)-y(i-1))./(2*h);
end

figure(1)
plot(t,y)
figure(2)
plot(t,dy)
%plot(t(2:n-1), df(2:n-1))
%c )

for i =2:200
    d2y(i) = y(i+1) - 2*y(i) + y(i-1)/(h^2);
end

pt = [t(2) t(3) t(4)];
py = [d2y(2) d2y(3) d2y(4)];
p = polyfit(pt,py,2);
d2y(1) = polyval(p, t(1));
pt = [t(198) t(199) t(200)];
py = [d2y(198) d2y(199) d2y(200)];
p = polyfit(pt, py,2);
d2y(201) = polyval(p,t(201));
figure(3)
plot(t,d2y)


%% uppg 2a
a = 0;
b = 3;

n = 201;
f = @(x) x.^2;
h = (b-a)/(n-1);

syms i
TC = double(symsum((h./2)*(fret(i) + fret(i+1)),i,1,(n-1)))
%% uppg 2b
nVal = [11 101 201]

for j = 1:1:3
    n = nVal(j);
    h = (3-0)./(n-1);
    itot = 0.0;
    for i = 2:2:(n-1)
        TC = sum(((h./3)*(fret(i+1) + 4*fret(i) + fret(i-1))),2:2:n-1);
        itot = itot + TC;
    end
    itot
end

%% uppg 3a

r = linspace(0,10);
p = @(r) 0.084./(1+exp((2*r)-8));
y = p(r);
plot(r,y)

%% uppg 3b

format long
a = 0;
b = 10000;
n = b;
r = linspace(a,b,n);
h = (b-a)/(n-1);

Q = 0.0;
y = @(r) (0.084./(1+exp((2*r)-8))).*(r.^2);

for i = 2:2:(n-1)
    TC = (h/3 * (y(i+1) + 4*y(i) + y(i-1)));
    Q = Q + TC;
end
Q = (4*pi*Q)

%% uppg 3c

Z = 10;


%% uppg4

a = 0;
b = 1;
c = 0;
d = 2;

m = (10);
n = (20);
f = @(x,y) x.^2 .*cos(y);

sum = sumComp(f,a,b,c,d,n,m)
sumquad = dblquad(f,a,b,c,d)

%% uppg5

% x : vektro som inneh??ller x-koordinater f??r alla punkter
% y : samma fast f??r y
% T : n x 3 matris, n = antal trianglar
% f(Ci) : central punkten

x = cos(0:0.3:2*pi);
y = sin(0:0.3:2*pi);

x = [0.01*x' ; 0.05*x' ; 1*x' ; 1.5*x' ; 2*x'];
y = [0.01*y' ; 0.05*y' ; 1*y' ; 1.5*y' ; 2*y'];

T = delaunay(x,y);
trimesh(T,x,y);
f = @(x,y) exp(-x.^2 -y.^2);
I = triangleSum(f,T,x,y)
disp('exakt v??rde'), disp(2.*(0.5-exp(-4)/2)*pi)

%% DEL2 uppg1a
x = linspace(0,10);
fdz = @(x) log(((x.^2)/2) + exp(-2));
figure(1)
y = fdz(x);
plot(x,y)
title("analytical")

dz = @(t,z) t*exp(-z);
[t,y] = ode45(dz,[0 10],-2);
figure(2)
plot(t,y)
title("ode45")

%% uppg 1b

%y'' + y = sin(x), y(0) = 0, y'(0) = 0

dy = @(x,y) [y(2); -y(1) + sin(x)];
[t,y] = ode45(dy,[0 10],[0 0]);
plot(t,y)
legend('-y(1) + sin(x)','y(2)')

%% uppg 2

[t,y] = ode45(@uppg2Func,[0 20], [10 0 15 0]');
plot(t,y(:,1),t,y(:,2),t,y(:,3),t,y(:,4))
legend('y1','dy1','y2','dy2')


%% uppg 3

m = 100;
v0 = 20;
k = 40;
g = 9.81;
time = [0 10];

dv = @(t,v) g - (k/m).*(v.^2);
[t,y] = ode45(dv,time,v0);
plot(t,y)

%% uppg 4

%Numerical a)
L = 1;
g = 9.81;

dy = @(x,y) [y(2); -(g/L).*sin(y(1))];
[t,y] = ode45(dy,[0 10],[pi/4 0]);
plot(t,y)
hold on

%Analytical b)
t = linspace(0,10);
fdy = @(t) (pi/4) .* cos((sqrt(g)*t)./(sqrt(L)));
y = fdy(t);
plot(t,y)

legend ('Numerical', 'Derrivative', 'Analytical')

%% uppg5a

c = 0.05;
g = 9.81;

dz = @(t,z) [z(2)
    -c*sqrt(z(2).^2 + z(4).^2) * z(2)
    z(4)
    -c*sqrt(z(2).^2 + z(4).^2)*z(4)-g];

opt = odeset('Event', @eventfun);
[t,z] = ode45(dz, [0 2], [0 10 0 10], opt);

plot(z(:,1),z(:,3))
grid on

%% uppg 5b

c = 0.25;
g = 9.81;
v = 10;

dz = @(t,z) [z(2)
    -c*sqrt(z(2).^2 + z(4).^2) * z(2)
    z(4)
    -c*sqrt(z(2).^2 + z(4).^2)*z(4)-g];

maxLen = -100000;
angle = 0;
for i = 1:89
    opt = odeset('Event', @eventfun);
    [t,z] = ode45(dz, [0 2], [0 v*cosd(i) 0 v*sind(i)], opt);
    
    if z(end,1) > maxLen
        maxLen = z(end,1);
        angle = i;
    end
end
[t,z] = ode45(dz, [0 2], [0 v*cosd(angle) 0 v*sind(angle)], opt);
plot(z(:,1),z(:,3))
grid on
angle
maxLen

%% uppg 5c

c = 0.25;
g = 9.81;
%v = 100;

dz = @(t,z) [z(2)
    -c*sqrt(z(2).^2 + z(4).^2) * z(2)
    z(4)
    -c*sqrt(z(2).^2 + z(4).^2)*z(4)-g];

maxLen = -100000;
angle = 0:10:10;
for v = 10:10:100
    for i = 1:89
        opt = odeset('Event', @eventfun);
        [t,z] = ode45(dz, [0 2], [0 v*cosd(i) 0 v*sind(i)], opt);
        
        if z(end,1) > maxLen
            maxx = t;
            maxy = z;
            maxLen = z(end,1);
            angle(v/10) = i;
        end
    end
    [t,z] = ode45(dz, [0 2], [0 v*cosd(angle(v/10)) 0 v*sind(angle(v/10))], opt);
    plot(z(:,1),z(:,3))
    hold on
    grid on
end
    ylim ([0 4])

legend("angle :" + num2str(angle(1)),"angle :" + num2str(angle(2)),"angle :" + num2str(angle(3)),"angle :" + num2str(angle(4)),"angle :" + num2str(angle(5)),"angle :" + num2str(angle(6)),"angle :" + num2str(angle(7)),"angle :" + num2str(angle(8)),"angle :" + num2str(angle(9)),"angle :" + num2str(angle(10)))
maxLen

% Ja vinkeln beroro p??

%% DEL3 uppg 1a

GM = 1;
dz = @(t,z) [z(2)
    -GM * (z(1) / (z(1).^2 + z(3).^2).^(3/2))
    z(4)
    -GM * (z(3) / (z(1).^2 + z(3).^2).^(3/2))];

for e = [0 0.5 0.9]
[t,z] = ode45(dz, [0 2*pi], [(1-e) 0 0 (((e+1) ./ (1-e)).^(1/2))]);   
plot(z(:,1),z(:,3))
hold on
grid on
end

legend('e = 0', 'e = 0.5', 'e = 0.9')

%% uppg 1b
axis equal
m = 1;
GM = 1;
dz = @(t,z) [z(2)
    -GM * (z(1) / (z(1).^2 + z(3).^2).^(3/2))
    z(4)
    -GM * (z(3) / (z(1).^2 + z(3).^2).^(3/2))];

e = [0 0.5 0.9];
tit = ["e = 0.0" "e = 0.5" "e = 0.9"];

figure('NumberTitle', 'off', 'Name', 'Kinetic energy');
for i = 1:1:3
[t,z] = ode45(dz, [0 2*pi], [(1-e(i)) 0 0 (((e(i)+1) ./ (1-e(i))).^(1/2))]');   
Ekin = (0.5) .* (m*(z(:,2).^2 + z(:,4).^2));
subplot(3,1,i);
plot(t, Ekin)
title(tit(i));
hold on
grid on
end

figure('NumberTitle', 'off', 'Name', 'Potential energy');
for i = 1:1:3
[t,z] = ode45(dz, [0 2*pi], [(1-e(i)) 0 0 (((e(i)+1) ./ (1-e(i))).^(1/2))]');   
Epot = (-GM * m) .* (1./sqrt((z(:,1).^2 + z(:,3).^2)));
subplot(3,1,i);
plot(t, Epot)
title(tit(i));
hold on
grid on
end

figure('NumberTitle', 'off', 'Name', 'Total energy');
for i = 1:1:3
[t,z] = ode45(dz, [0 2*pi], [(1-e(i)) 0 0 (((e(i)+1) ./ (1-e(i))).^(1/2))]');   
Etot = (((z(:,2).^2 + z(:,4).^2) / 2) - (1 ./ sqrt((z(:,1).^2 + z(:,3).^2))));
subplot(3,1,i);
plot(t, Etot)
title(tit(i));
hold on
grid on
end

%% uppg 1c


GM = 1;
dz = @(t,z) [z(2)
    -GM * (z(1) / (z(1).^2 + z(3).^2).^(3/2))
    z(4)
    -GM * (z(3) / (z(1).^2 + z(3).^2).^(3/2))];

e = 0.5;

for a = 1:0.05:1.5
[t,z] = ode45(dz, [0 5*pi], [(1-e) 0 0 (((e+a) ./ (1-e)).^(1/2))]);   
plot(z(:,1),z(:,3))
xlim([-3 1])
hold on
grid on
end

legend('e = 0', 'e = 0.5', 'e = 0.9')



%% uppg 2a
%w0 = sqrt(k/m)

%% uppg 2b
m = 1;
k = 1;
F0 = 1;
w = sqrt(k/m)

dz = @(t,z) [z(2); -(k*z(1))/m + (F0 * cos(w*t))/m];
[t,z] = ode45(dz, [0 50], [0 0]);   
plot(t,z(:,1))
grid on

%% uppg 2c

m = 1;
k = 1;
F0 = 1;
c = 0.1;
w = sqrt(k/m)

dz = @(t,z) [z(2); ((c*z(2))/m) - (k*z(1))/m + (F0 * cos(w*t))/m];
[t,z] = ode45(dz, [0 50], [0 0]);   
plot(t,z(:,1))
grid on

%% uppg 2d

m = 1;
k = 1;
F0 = 1;
%c = 1;

for c = 0:1:10
t = linspace(1,2);
opt = odeset('initialStep',0.0033,'MaxStep',0.0033);
Fext = (0.5*(((abs(t-1)) / (t-1)) - ((abs(t-2)) / (t-2))))/m;
dz = @(t,z) [z(2); -((c*z(2))/m) - ((k*z(1))/m) + Fext];
[t,z] = ode45(dz, [0 10], [0 0]);
subplot(5,3,c + 1);
plot(t,z(:,1))
title(['c =  ',num2str(c)])
grid on
end


%% uppg 3a

format long
lambda1 = log(2) / 5.01;
lambda2 = log(2) / 138.38;

dz = @(t,z) [-lambda1*z(1)
    lambda1*z(1) - (lambda2 * z(2))];
[t,z] = ode45(dz, [0 100], [1 0]);
plot(t,z)
grid on

%b
%74.34 sekunder blir den max
%% uppg 4
k = 1;

gamma = fzero(@shooting,25);
dT = @(x,T) [T(2) ;  -(x-3+5*sin(pi*x))/k];
[x,T] = ode45(dT, [0 10],[15 ; gamma]);
plot(x,T(:,1))

%% uppg 5

a = input('a bound ');
b = input('b bound ');
y0 = a;
n = input('number of points ');
ab = (b - a);
h = ab / n;
f = @(x,y) -y + sin(x);
t0 = a;

[x,T] = ode45(f, [a b],1,n);
plot(x,T(:,1))

Eulers(f,a,b,y0,h,n,t0)

