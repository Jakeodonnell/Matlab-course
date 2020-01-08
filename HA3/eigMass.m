function eigMass = eigMass(masspoints,nodes)
massPoints = masspoints;
node = nodes;
time = 3;
twoVector = ones(1,massPoints)*-2;
oneVector = ones(1,(massPoints-1))*1;

twos_diagnal = diag(twoVector);
ones_above_diagnal = diag(oneVector,1);
ones_under_diagnal = diag(oneVector,-1);

complete_matrix = twos_diagnal + ones_above_diagnal + ones_under_diagnal;
[C,D] = eig(complete_matrix); % Eigenvectors and eigenvalues from complete matrix

ending = massPoints - 1;
l = [0:ending]'; %  Positions of the end points and the mass points, column vector

figure(1)
k = sqrt(-D(node,node));
for t = 0:0.1:(time*2)*pi % Duration of how long animation should go on
    u = sin(k*t)*[C(:,node)];% u as a column vector,
    % add positions for the end point
    % u = sin(k*t)*[0 ; C(:,1) ; 0]; % u as a column vector,
    plot(l,u,l,u, 'o')
    ylim([-1 1])
    pause (0.1)
end
end

