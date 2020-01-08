function deltaT = shooting(gamma)
gamma
k = 1;
dT = @(x,T) [T(2) ;  -(x-3+5*sin(pi*x))/k];
[x,T] = ode45(dT, [0 10],[15 ; gamma]);
deltaT = 15 - T(length(x),1)
end

