function [M1 M2] = dens()
%DENS Summary of this function goes here
%   Detailed explanation goes here

%a--------------
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
Mass = (Nsphere / N)*8*R.^3* fsum/Nsphere


%b--------------
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
            fsum = fsum + 14000 * (x.^2 + y.^2);
        elseif r2 < 3460000
            Nsphere = Nsphere + 1;
            fsum = fsum + 11000* (x.^2 + y.^2);
        elseif r2 < 5630000
            Nsphere = Nsphere + 1;
            fsum = fsum + 4850* (x.^2 + y.^2);
        elseif r2 < 6340000
            Nsphere = Nsphere + 1;
            fsum = fsum + 3800* (x.^2 + y.^2);
        elseif r2 < 6370000
            Nsphere = Nsphere + 1;
            fsum = fsum + 2600 * (x.^2 + y.^2);
        end
    end
end

realInertia = (Nsphere / N)*8*R.^3 * fsum/Nsphere
calcmass_inertia = (2/5) * Mass*(R^2)
constant_inertia = (2/5) * (5514 * (4*pi)/3 * R^3) * R.^2

%Mass_with_cesity_constant = 5514 * (4*pi/3) * R^3
%Inertia_with_cesity_constant = 2/5 * Mass_with_cesity_constant * R^2
