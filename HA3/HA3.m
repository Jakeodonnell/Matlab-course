%% 2
eigMass(3,3)
%% 3
eigMass(10,1)
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

