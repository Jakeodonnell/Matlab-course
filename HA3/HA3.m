%% 2
eigMass(3,1)
%% 3
eigMass(10,2)
%% uppg 6
% steg (b-a) / (n-1)
% y = [y2, y3... yn-1]';
% b = (f2 + (y(a)/ delta(x^2), f3,f4,..., f(n-1) + y(b)/deltax^2  ))

a = 0;
b = 10;

Ta = 15;
Tb = 15;

k = 1;
n = 500;

x = linspace(a,b,n)'; %kolonvector
f = x - 3 + 5*sin(pi*x); % f1,f2,f3,f4...fn

dx = (b-a) / (n-1);
% skapa 2
L = zeros(n-2);
for i = 1:n-2
    L(i,i) = 2 * (1/dx^2);
end

for i = 1:n-3
    L(i,i+1) = -1/(dx^2);
    L(i+1,i) = -1/(dx^2);
end

B = [Ta/dx^2; zeros(n-4,1); Tb/dx^2] + (1/k) * f(2:n-1);

T = L\B;
T = [Ta;T;Tb];

figure(2)
plot(x,T),xlabel('x'),ylabel('y')
grid on

%% DEL 2 uppg 1
a = 0; % Start of interval
b = 10; % End of interval
N = 100000; % Number of random numbers
Random_numbers = [10^1 10^2 10^3 10^4 10^5 10^6 10^7 10^8];
results = [];

 for c = 1:1:8
        N = Random_numbers(c);
        x = a+(b-a)*rand(1,N); % A vector x is constructed with 100 000 random
        % numbers on the interval [0, 2pi]
        y = x.*exp(-x); % A vector y is constructed with 100 000 random
        % function values on the interval [0, 2pi]
        Iest = (b-a)/N*sum(y); % Computes an estimate of the integral
        results(c) = Iest;
 end
Random_numbers
results

%% uppg 2

N = 100000;
R = 6370000;
fsum = 0;
Nsphere = 0;

for i=1:N
    x = -R+(2*R*rand);                % x is assigned a random number between -R and R
    y = -R+(2*R*rand);                % y is assigned a random number between -R and R
    z = -R+(2*R*rand);                % z is assigned a random number between -R and R
    
    r2 = sqrt(x.^2 + y.^2 + z.^2);      % Distance squared of the point (x,y,z) from origo 
    
    if r2 < R
        Nsphere = Nsphere + 1;
        if r2 < 1200000
            Nsphere = Nsphere + 1;
            fsum = fsum + 14000;
        elseif r2 < 3460000
             Nsphere = Nsphere + 1;
             fsum = fsum + 11000;
        elseif r2 < 5630000
             Nsphere = Nsphere + 1;
             fsum = fsum + 4850;
        elseif r2 < 6340000
             Nsphere = Nsphere + 1;
             fsum = fsum + 3800;
        elseif r2 < 6370000
             Nsphere = Nsphere + 1;
             fsum = fsum + 2600;
        end
    end
end
size = (Nsphere / N)*8*R.^3* fsum/Nsphere