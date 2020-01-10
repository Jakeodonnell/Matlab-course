figure(1)
clf
E0 = 0.7;                   % Initial energy
npart = 1000;                % Number of particles to be simulated
absprob = 0.2;              % Probability for interaction through absorbtion
arrWithEdep = [];
% Construction of circular target
%r = [0.1 0.2 0.3 0.4 0.5];                    % Target radius
for a = 1:1:5
    r = a/10;
    center = 0.5;               % Center of target on y-axis
    th = linspace(0,2*pi);
    xx = r*cos(th);
    yy = center+r*sin(th);
    plot(xx,yy,'r');
    hold on
    
    Edep = 0;                   % Initialize deposited energy
    for j=1:npart               % Loop over number of particles
        clear xtrack             % Clear track in x-direction
        clear ytrack             % Clear track in y-direction
        
        i = 1;                   % Initialize interaction counter i
        k = 0;
        
        % Position of particle gun and initial particle energy
        E = E0;                  % Initial energy of particles
        xtrack(i) = 0;           % Particles are fired from x=0 and
        ytrack(i) = 0;           % y=0
        
        i = i + 1;
        
        d = (center-r) + E*rand;    % Deduce range of particles before interaction
        if d < (center+r)           % If particle interacts within target
            xtrack(i) = xtrack(i-1); % Update x-position of particle track
            ytrack(i) = ytrack(i-1) + d; % Update y-position of particle track
            thetaold = 0;            % Store direction of particle
            i = i + 1;
        else                        % Particle does not interact within target
            xtrack(i) = xtrack(i-1); % Particle is transported far outside the target
            ytrack(i) = 1000*d;
            k = 1;                   % End track
        end
        
        while k<1                   % While particle is inside target and not yet absorbed
            pint = rand;
            if pint < absprob        % If absorbition
                Edep = Edep + E;      % Increment deposited energy
                k = 1;                % End particle track
            else                     % If elastic scattering
                theta = -pi+2*pi*rand + thetaold;  % Determine deflection angle
                E = 0.5*E;            % Update energy
                Edep = Edep + E;      % Increment deposited energy
                d = E*rand;           % Range before interaction
                x = d*sin(theta);     % Change in x-position due to scattering
                y = d*cos(theta);     % Change in y-position due to scattering
                thetaold = theta;     % Store direction of particle
                xtrack(i) = xtrack(i-1) + x;  % Update x-position of track
                ytrack(i) = ytrack(i-1) + y;  % Update y-position of track
                
                % If new position outside target the particle is transported
                % away from the target and track ends
                if (xtrack(i).^2+(ytrack(i)-center).^2) > r.^2
                    xtrack(i) = 1000*x;
                    ytrack(i) = 1000*y;
                    k = 1;                  % End particle track
                end
                i = i + 1;
            end
        end
        subplot(3,2,(r*10))
        plot(xtrack,ytrack,'--o')        % Plot particle track
        hold on
    end
    xlabel('x');
    ylabel('y');
    axis([-0.6 0.6 0 1])
    Edep = Edep/npart                   % Compute average deposited energy
    arrWithEdep(r*10) = Edep;
    title(Edep)
end
arrEdep = arrWithEdep
averageOfAll = median(arrWithEdep)
%% b
figure(2)

clf
E0 = 0.7;                   % Initial energy
npart = 1000;                % Number of particles to be simulated
% Construction of circular target
% Target radius
%absprob = 0.2;              % Probability for interaction through absorbtion
arrWithAbs = [];
AbsValues = [0.1 0.3 0.5 0.7 0.9];
ind = 0;
for z = 1:1:5
    absprob = AbsValues(z)
    r = 0.3;
    center = 0.5;               % Center of target on y-axis
    th = linspace(0,2*pi);
    xx = r*cos(th);
    yy = center+r*sin(th);
    plot(xx,yy,'r');
    hold on
    
    Edep = 0;                   % Initialize deposited energy
    for j=1:npart               % Loop over number of particles
        clear xtrack             % Clear track in x-direction
        clear ytrack             % Clear track in y-direction
        
        i = 1;                   % Initialize interaction counter i
        k = 0;
        
        % Position of particle gun and initial particle energy
        E = E0;                  % Initial energy of particles
        xtrack(i) = 0;           % Particles are fired from x=0 and
        ytrack(i) = 0;           % y=0
        
        i = i + 1;
        
        d = (center-r) + E*rand;    % Deduce range of particles before interaction
        if d < (center+r)           % If particle interacts within target
            xtrack(i) = xtrack(i-1); % Update x-position of particle track
            ytrack(i) = ytrack(i-1) + d; % Update y-position of particle track
            thetaold = 0;            % Store direction of particle
            i = i + 1;
        else                        % Particle does not interact within target
            xtrack(i) = xtrack(i-1); % Particle is transported far outside the target
            ytrack(i) = 1000*d;
            k = 1;                   % End track
        end
        
        while k<1                   % While particle is inside target and not yet absorbed
            pint = rand;
            if pint < absprob        % If absorbition
                Edep = Edep + E;      % Increment deposited energy
                k = 1;                % End particle track
            else                     % If elastic scattering
                theta = -pi+2*pi*rand + thetaold;  % Determine deflection angle
                E = 0.5*E;            % Update energy
                Edep = Edep + E;      % Increment deposited energy
                d = E*rand;           % Range before interaction
                x = d*sin(theta);     % Change in x-position due to scattering
                y = d*cos(theta);     % Change in y-position due to scattering
                thetaold = theta;     % Store direction of particle
                xtrack(i) = xtrack(i-1) + x;  % Update x-position of track
                ytrack(i) = ytrack(i-1) + y;  % Update y-position of track
                
                % If new position outside target the particle is transported
                % away from the target and track ends
                if (xtrack(i).^2+(ytrack(i)-center).^2) > r.^2
                    xtrack(i) = 1000*x;
                    ytrack(i) = 1000*y;
                    k = 1;                  % End particle track
                end
                i = i + 1;
            end
        end
        subplot(3,2,z)
        plot(xtrack,ytrack,'--o')        % Plot particle track
        hold on
    end
    xlabel('x');
    ylabel('y');
    axis([-0.6 0.6 0 1])
    Edep = Edep/npart;                   % Compute average deposited energy
    arrWithAbs(z) = Edep;
end
arrAbs = arrWithAbs
averageOfAll = median(arrWithAbs)

%% c
figure(3)

clf
E0 = 0.7;                   % Initial energy
npart = 1000;                % Number of particles to be simulated
absprob = 0.2;
r = 0.3;
center = 0.5;               % Center of target on y-axis
th = linspace(0,2*pi);
xx = r*cos(th);
yy = center+r*sin(th);
plot(xx,yy,'r');
hold on
Edep = 0;                   % Initialize deposited energy
for j=1:npart               % Loop over number of particles
    clear xtrack             % Clear track in x-direction
    clear ytrack             % Clear track in y-direction
    
    i = 1;                   % Initialize interaction counter i
    k = 0;
    
    % Position of particle gun and initial particle energy
    E = E0;                  % Initial energy of particles
    xtrack(i) = 0;           % Particles are fired from x=0 and
    ytrack(i) = 0;           % y=0
    
    i = i + 1;
    
    d = (center-r) + E*rand;    % Deduce range of particles before interaction
    if d < (center+r)           % If particle interacts within target
        xtrack(i) = xtrack(i-1); % Update x-position of particle track
        ytrack(i) = ytrack(i-1) + d; % Update y-position of particle track
        thetaold = 0;            % Store direction of particle
        i = i + 1;
    else                        % Particle does not interact within target
        xtrack(i) = xtrack(i-1); % Particle is transported far outside the target
        ytrack(i) = 1000*d;
        k = 1;                    % End track
    end
    
    while k<1                     % While particle is inside target and not yet absorbed
        pint = rand;
        if pint < absprob         % If absorbition
            Edep = Edep + E;      % Increment deposited energy
            k = 1;                % End particle track
        else                      % If elastic scattering
            phi = -pi+2*pi*rand;
            theta = phi/2;        % Determine deflection angle
            E = 0.5*E*(1-cos(phi));                % Update energy
            Edep = Edep + E;      % Increment deposited energy
            d = E*rand;           % Range before interaction
            x = d*sin(theta);     % Change in x-position due to scattering
            y = d*cos(theta);     % Change in y-position due to scattering
            thetaold = theta;     % Store direction of particle
            xtrack(i) = xtrack(i-1) + x;  % Update x-position of track
            ytrack(i) = ytrack(i-1) + y;  % Update y-position of track
            
            % If new position outside target the particle is transported
            % away from the target and track ends
            if (xtrack(i).^2+(ytrack(i)-center).^2) > r.^2
                xtrack(i) = 1000*x;
                ytrack(i) = 1000*y;
                k = 1;                  % End particle track
            end
            i = i + 1;
        end
    end
    plot(xtrack,ytrack,'--o')        % Plot particle track
    hold on
end
xlabel('x');
ylabel('y');
axis([-0.6 0.6 0 1])
Edep = Edep/npart                   % Compute average deposited energy

%% d
figure(4)

clf
E0 = 0.7;                   % Initial energy
npart = 1000;                % Number of particles to be simulated
absprob = 0.1;
r = 0.3;
center = 0.5;               % Center of target on y-axis
th = linspace(0,2*pi);
xx = r*cos(th);
yy = center+r*sin(th);
plot(xx,yy,'r');
hold on
o = 1;
Edep = 0;                   % Initialize deposited energy

print = false;
counter = 0;

for j=1:npart               % Loop over number of particles
    clear xtrack             % Clear track in x-direction
    clear ytrack             % Clear track in y-direction
    
    i = 1;                   % Initialize interaction counter i
    k = 0;
    
    % Position of particle gun and initial particle energy
    E = E0;                  % Initial energy of particles
    xtrack(i) = 0;           % Particles are fired from x=0 and
    ytrack(i) = 0;           % y=0
    
    i = i + 1;
    
    d = (center-r) + E*rand;    % Deduce range of particles before interaction
    if d < (center+r)           % If particle interacts within target
        xtrack(i) = xtrack(i-1); % Update x-position of particle track
        ytrack(i) = ytrack(i-1) + d; % Update y-position of particle track
        thetaold = 0;            % Store direction of particle
        i = i + 1;
    else                        % Particle does not interact within target
        xtrack(i) = xtrack(i-1); % Particle is transported far outside the target
        ytrack(i) = 1000*d;
        k = 1;                   % End track
    end
    
    while k<1                   % While particle is inside target and not yet absorbed
        pint = rand;
        if pint < absprob        % If absorbition
            Edep = Edep + E;      % Increment deposited energy
            k = 1;                % End particle track
        else                     % If elastic scattering
            ra = uint8(rand*10);
            if (mod(ra,2) == 0)
                print = true;
            end
            phi = -pi+2*pi*rand;
            theta = phi/2;  % Determine deflection angle
            E = 0.5*E*(1-cos(phi));            % Update energy
            Edep = Edep + E;      % Increment deposited energy
            d = E*rand;           % Range before interaction
            x = d*sin(theta);     % Change in x-position due to scattering
            y = d*cos(theta);     % Change in y-position due to scattering
            thetaold = theta;     % Store direction of particle
            xtrack(i) = xtrack(i-1) + x;  % Update x-position of track
            ytrack(i) = ytrack(i-1) + y;  % Update y-position of track
            
            % If new position outside target the particle is transported
            % away from the target and track ends
            if (xtrack(i).^2+(ytrack(i)-center).^2) > r.^2
                xtrack(i) = 1000*x;
                ytrack(i) = 1000*y;
                k = 1;                  % End particle track
            end
            i = i + 1;
        end
    end
    if(print && counter < 10)
        plot(xtrack,ytrack,'--o')        % Plot particle track
        hold on
        print = false;
        counter = counter + 1;
    end
end
xlabel('x');
ylabel('y');
axis([-0.6 0.6 0 1])
Edep = (Edep/npart);
title(Edep)
