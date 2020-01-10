function t = eigMass(masspoints)
massPoints = masspoints;
node = [1,2,3,4,5,6,7,8,9,10];
time = 2;
twoVector = ones(1,massPoints)*-2;
oneVector = ones(1,(massPoints-1))*1;

twos_diagnal = diag(twoVector);
ones_above_diagnal = diag(oneVector,1);
ones_under_diagnal = diag(oneVector,-1);

complete_matrix = twos_diagnal + ones_above_diagnal + ones_under_diagnal;
[C,D] = eig(complete_matrix); % Eigenvectors and eigenvalues from complete matrix

ending = massPoints +1;
l = [0:ending]'; %  Positions of the end points and the mass points, column vector

figure(1)
for i = 1:1:masspoints
    clf
    k = sqrt(-D(node(i),node(i)));
    for t = 0:0.1:(time)*pi % Duration of how long animation should go on
        u = sin(k*t)*[0; C(:,node(i)); 0];% u as a column vector,
        % add positions for the end point
        % u = sin(k*t)*[0 ; C(:,1) ; 0]; % u as a column vector,
        plot(l,u,l,u, 'o')
        ylim([-1 1])
        pause (0.01)
    end
end
end

