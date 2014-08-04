%
% Finite difference discretisation of the equation
%          -epsilon y'' + y' = 2
% using central difference
%
clear all; close all; clc
epsilon = 0.01;
xsol  = (0:0.01:1);
% Analytical solution:
ysol = 2*xsol + 2*(exp(-1/epsilon) - exp(-(1-xsol)/epsilon))/(1 - exp(-1/epsilon));

scrsz = get(0,'ScreenSize');
figHandle1 = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/1]);

lim = 205;
ns = [4:10 15:10:lim];
lim_ns = length(ns);
for i = 1:1:lim_ns
    n=ns(i);
    h = 1/n;
    % N: number of internal gridpoints (number of unknowns)
    N = n-1;
    
    % Subdiagonal elements:
    c1= -epsilon/h^2 - 1/(2*h);
    % Main diagonal elements:
    c2= 2*epsilon/h^2;
    % Super diagonal elements:
    c3= -epsilon/h^2 + 1/(2*h);
    % Right-hand side vector:
    b = 2*ones(N,1);
    % Use matlab sparse storage (only nonzeros are stored): 
    A = gallery('tridiag',N,c1,c2,c3);
    % Compute solution:
    w = A\b;
    % Add the boundary conditions for plotting:
    y = zeros(N+2,1);
    y(2:N+1) = w';
    x = (0:h:1);
    % Plot (in same figure as analytical solution:

    figure(figHandle1);
        

    subplot(2,2,[1 2 ])
    plot(xsol,ysol)
    hold on;
    plot(x,y,'-*r');
    legend('Solution','Approximation','Location','BestOutside')
    title(sprintf('Finite diference discretisation Conv-Diff Eq. (central diff) #points = %d',n))
    axis([-0.01 1.01 -3 6])
    
    ea = eig(full(A));
    subplot(2,2,[3 4])
    circle([c2,0],abs(c1)+abs(c3),100,'k'); hold on; axis equal;
    plot(real(ea),imag(ea),'*')
    line([0 0], [abs(c1)+abs(c3) -(abs(c1)+abs(c3))]);
    legend('Gershgorin circle','Eigenvalues','Location','BestOutside')
    title('Eigenvalues location and Gershgorin circle')
     
    pause(0.3);
    if i< lim_ns
        clf(figHandle1)
        end
    
end