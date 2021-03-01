%
% Khaled AbouShaban 
% Elec4700 Assignment 2 
% Q1
%

clear all
close all
set(0, 'DefaultFigureWindowStyle', 'docked')


W = 2;
L = 3;
V0 = 1;

dx = 0.2;
dy = 0.2;
nx = 50;
ny = L/W *50;


%
% implementing using matrix
%
C1 = -2*(1/dx^2 + 1/dy^2);
C2 = 1/(dx^2);
C3 = 1/(dy^2);
G = zeros(nx*ny,nx*ny);
Y = zeros(nx*ny,1);

for x=2:(nx-1)
    for y=2:(ny-1)
        i = (y-1).*nx + x;
        G(i,i) = C1;
        G(i,((y-1).*nx + x-1)) = C2;
        G(i,((y-1).*nx + x+1)) = C2;
        G(i,((y-2).*nx + x)) = C3;
        G(i,((y).*nx + x)) = C3;
    end
end


for y=1:ny
    i = ((y-1).*nx + 1);
    G(i,i) = 1;
    
    Y(i) = V0;
    
    i = ((y-1).*nx + nx);
    G(i,i) = 1;
end


for x=2:(nx-1)
    i = (1-1).*nx + x;
    G(i,i) = 1;
    G(i,(2-1).*nx + x) = -1;
    
    i = ((ny-1).*nx + x);
    G(i,i) = 1;
    G(i,((ny-2).*nx + x)) = -1;
end
V = G\Y;
V = reshape(V,[],ny)';

figure();
surf(linspace(0,L,nx),linspace(0,W,ny),V);
xlabel('x')
ylabel('y')
title(sprintf('2-D plot of V(x) = %.2f ', dx))
set(gca, 'View', [45 45])



analy = zeros(ny, nx);
x1 = repmat(linspace(-L/2,L/2,nx),ny,1);
y1 = repmat(linspace(0,W,ny),nx,1)';
itr = 150;
avgError = zeros(itr,1);


for i=1:itr
    n = 2*i - 1;
    analy = analy + 1./n.*cosh(n.*pi.*x1./W) ...
        ./cosh(n.*pi.*(L./2)./W).*sin(n.*pi.*y1./W);

    avgError(i) = mean(mean(abs(analy.*4.*V0./pi - V)));
end

analy = analy.*4.*V0./pi;

figure();
surf(linspace(0,L,nx),linspace(0,W,ny),analy);
xlabel('x')
ylabel('y')
title(sprintf('Analytical Graph with %d iterations ', itr))