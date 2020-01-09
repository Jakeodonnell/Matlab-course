% Program by Ian Cooper, School of Physics, University of Sydney
% for solving time independent Schrodinger Equation for a potential well
% Modified and adapted by Per Jonsson, 2016.

clear all
close all
clc
num = 3001; % Number of grid points (odd number)
clf
set(gcf, 'Position',  [2000, 50, 800, 800])
% Constants ------------------------------------------------------------
hbar = 1.055e-34; % J.s
me = 9.109e-31; % kg electron mass
Ese = 1.6e-19; % Energy scaling factor
Lse = 1e-9; % Length scaling factor
Cse = -hbar^2/(2*me) / (Lse^2*Ese); % Schrodinger Eq constant

% Potential well parameters, well centered around x = 0 -----------------
% Width 0.1 nmm depth - 400 eV

U = zeros(num,1);
U2 = zeros(num,2);

U_matrix = zeros(num-2);
xMin = -0.2; % default = -0.25 nm
xMax = 0.2; % default = +0.25 nm
x1 = 0.15;
x2 = 0.05; % 1/2 well width: default = 0.05 nm
U1 = -400; % Depth of well (eV): default = -400 eV
U2 = -100;
x = linspace(xMin,xMax,num);
for cn = 1 : num
    if abs(x(cn)) <= x1
        U(cn) = U2;
    end
    if abs(x(cn)) <= x1 && abs(x(cn)) >= x2
        U(cn) = U1;
    end;
end;

% Plot potential well --------------------------------------------
figure(1);
subplot(3,4,1);
plot(x,U,'LineWidth',3);
grid on
axis([xMin-eps xMax min(U)-50 max(U)+50]);
title('Square well','FontSize',12)
xlabel('x (nm)','FontSize',12);
ylabel('energy (eV)','FontSize',12);


% Make potential energy matrix ----------------------------------
dx = (x(2)-x(1));
dx2 = dx^2;
for cn =1:(num-2)
    U_matrix(cn,cn) = U(cn+1); % diagonal matrix
end;

% Make Second Derivative Matrix ------------------------------------------
off = ones(num-3,1);
SD_matrix = (-2*eye(num-2) + diag(off,1) + diag(off,-1))/dx2;

% Make KE Matrix
K_matrix = Cse*SD_matrix;

% Make Hamiltonian Matrix
H_matrix = K_matrix + U_matrix;

% Find Eignevalues E_n and Eigenfunctions psi_N ---------------------------

[e_funct e_values] = eig(H_matrix);

% All Eigenvalues 1, 2 , ... nbound where E < 0 (bounded states)
flag = 0;
n = 1;
while flag == 0
    E(n) = e_values(n,n);
    if E(n) > 0
        flag = 1;
    end
    n = n + 1;
end
nbound = n-2;

% Corresponding Eigenfunctions 1, 2, ... ,nbound
for cn = 1 : nbound
    psi(:,cn) = [0; e_funct(:,cn); 0]; % add zeros at the end points
    prob(:,cn) = psi(:,cn).* psi(:,cn); % probability
    if psi(5,cn) < 0
        psi(:,cn) = -psi(:,cn);
    end % curve starts positive
end % for

% Display eigenvalues in Command Window -------------------------------
disp(' ');
disp('================================================================ ');
disp(' ');
fprintf('No. bound states found = %0.0g \n',nbound);
disp(' ');
disp('Quantum State / Eigenvalues En (eV)');
for cn = 1 : nbound
    fprintf(' %0.0f ',cn);
    fprintf(' %0.5g \n',E(cn));
end
disp(' ')
disp(' ');

% Plot energy spectrum ------------------------------------------------
xs(1) = -x1;
xs(2) = x1;
%figure(2);
subplot(3,4,2);
plot(x,U,'b','LineWidth',2);
xlabel('position x (nm)','FontSize',12);
ylabel('energy U, E_n (eV)','FontSize',12);
title('Square well with energies','FontSize',12);
hold on
for cn = 1 : nbound
    ys(1) = E(cn);
    ys(2) = ys(1);
    plot(xs,ys,'r','LineWidth',2);
end %for
axis([xMin-eps xMax min(U)-50 max(U)+50]);

% Plots the eigenfunctions (wavefunctions)
for cn = 1:nbound
    %figure(cn+2)
    subplot(3,4,cn+2);
    y1 = psi(:,cn)/(max(psi(:,cn)-min(psi(:,cn))));
    y2 = prob(:,cn)/max(prob(:,cn));
    y3 = 1 + 2*U/(max(U) - min(U));
    plot(x,y1,'r','lineWidth',2)
    hold on
    plot(x,y2,'k','lineWidth',2)
    plot(x,y3,'b','lineWidth',2)
    axis off
    title_m = ['\psi n = ', num2str(cn)] ;
    title(title_m,'Fontsize',10);
end
