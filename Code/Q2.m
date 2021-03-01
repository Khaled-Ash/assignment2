%
% Elec4700 Q2
%

% Solving Current Flow with the Finite Diffrenece Method FDM in rectangluar
% region

clear all
close all
set(0, 'DefaultFigureWindowStyle', 'docked')

nx = 75;
ny = 50;
L = 20;
W = 10;
V1 = 1; 


SigmaConduct = 1;
SigmaInsulate = 10e-2;
CM = SigmaConduct*ones(nx, ny);
CM(1:W,(1:L)+ny/2-L/2) = SigmaInsulate;
CM((1:W)+nx-W,(1:L)+ny/2-L/2) = SigmaInsulate;

figure();
hold on;
surf(linspace(0,1.5,ny), linspace(0,1,nx), CM,'EdgeColor','none','LineStyle','none');
xlabel('x');
ylabel('y');
zlabel('Conduction (Mho)');
view([120 25])
title('Sigma(x,y)')

% Part A

V = NuS(nx, ny, CM, Inf, Inf, 0, V1);

figure();
hold on;
surf(linspace(0,1.5,ny), linspace(0,1,nx), V,'EdgeColor','none','LineStyle','none');
xlabel('x');
ylabel('y');
zlabel('Voltage (V)');
view([120 25])
colorbar
title('V(x,y)')

[Ex, Ey] = gradient(V);
Ex = -Ex;
Ey = -Ey;

figure();
quiver(linspace(0,1.5,ny), linspace(0,1,nx), Ex, Ey);
ylim([0 1]);
xlim([0 1.5]);
xlabel('x');
ylabel('y');
title('E(x,y)')

Jx = CM.*Ex;
Jy = CM.*Ey;
J = sqrt(Jx.^2 + Jy.^2);

figure();
hold on;
contourf(linspace(0,1.5,ny), linspace(0,1,nx), J,'EdgeColor','none','LineStyle','none');
quiver(linspace(0,1.5,ny), linspace(0,1,nx), Jx, Jy);
xlabel('x');
ylabel('y');
colorbar
title('J(x,y)')

% Part B

figure();
hold on;
grid on
range = 20:10:100;
I = [];
for x = range
    
 I = [I TotalCur(x, ny, V1, SigmaConduct, SigmaInsulate, W, L)];
    
end
plot(range, I);
ylabel('Current (A)');
xlabel('Mesh size');
title('Plot of Total Current vs Mesh Width ')

% Part C

figure();
range = 0:1:50;
I = [];
for w = range
    
    I = [I TotalCur(nx, ny, V1, SigmaConduct, SigmaInsulate, w, L)];
end
plot(range, I);
ylabel('Current (A)');
xlabel('Width');
title('Plot of Total current vs Box Width')
grid on 

% Part D

figure();
hold on;
grid on
range = logspace(-5,0, 50);
I = [];
for sigma = range
    
    I = [I TotalCur(nx, ny, V1, SigmaConduct, sigma, W, L)];

end
plot(range, I);
ylabel('Current (A)');
xlabel('Conduction (Mho)');
title('Plot of Total Current vs Box Conduction')

% Functions

function V = NuS(nx, ny, c, BLeft, BRight, BTop, BBottom)

    global C;
    g = sparse(nx*ny, ny*nx);
    b = zeros(1, nx*ny);
    for x=1:nx
        for y=1:ny
            n = y + (x - 1)*ny;
            nxm = y + (x - 2)*ny;
            nxp = y + (x)*ny;
            nym = (y-1) + (x - 1)*ny;
            nyp = (y+1) + (x - 1)*ny;

            if (x == 1 && y == 1)
                if (BLeft == Inf)
                    rxp = (c(x,y) + c(x+1,y))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;

                    g(n,n)   = -(rxp+ryp);
                    g(n,nxp) =  rxp;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BLeft;
                end
            elseif (x == 1 && y == ny)
                if (BLeft == Inf)
                    rxp = (c(x,y) + c(x+1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;

                    g(n,n)   = -(rxp+rym);
                    g(n,nxp) =  rxp;
                    g(n,nym) =  rym;
                else
                    g(n,n) = 1;
                    b(n) = BLeft;
                end
            elseif x == nx && y == 1 
                if (BRight == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;
                    g(n,n)   = -(rxm+ryp);
                    g(n,nxm) =  rxm;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BRight;
                end
            elseif x == nx && y == ny
                if (BRight == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;
                    g(n,n)   = -(rxm+rym);
                    g(n,nxm) =  rxm;
                    g(n,nym) =  rym;
                else
                    g(n,n) = 1;
                    b(n) = BRight;
                end
            elseif (x == 1)
                if (BLeft == Inf)
                    rxp = (c(x,y) + c(x+1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;

                    g(n,n)   = -(rxp+rym+ryp);
                    g(n,nxp) =  rxp;
                    g(n,nym) =  rym;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BLeft;
                end
                
            elseif x == nx 
                if (BRight == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;
                    g(n,n)   = -(rxm+rym+ryp);
                    g(n,nxm) =  rxm;
                    g(n,nym) =  rym;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BRight;
                end
               
            elseif y == 1
                if (BTop == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    rxp = (c(x,y) + c(x+1,y))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;
                    g(n,n) = -(rxm+rxp+ryp);
                    g(n,nxm) =  rxm;
                    g(n,nxp) =  rxp;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BTop;
                end
                
            elseif y == ny 
                if (BBottom == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    rxp = (c(x,y) + c(x+1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;
                    g(n,n) = -(rxm+rxp+rym);
                    g(n,nxm) =  rxm;
                    g(n,nxp) =  rxp;
                    g(n,nym) =  rym;
                else
                    g(n,n) = 1;
                    b(n) = BBottom;
                end
            else
                rxm = (c(x,y) + c(x-1,y))/2;
                rxp = (c(x,y) + c(x+1,y))/2;
                rym = (c(x,y) + c(x,y-1))/2;
                ryp = (c(x,y) + c(x,y+1))/2;
                
                g(n,n)   = -(rxm+rxp+rym+ryp);
                g(n,nxm) =  rxm;
                g(n,nxp) =  rxp;
                g(n,nym) =  rym;
                g(n,nyp) =  ryp;
            end
        end
    end
    
    V_temp = g\b';
    
    V = zeros(nx,ny,1);
    for x=1:nx
        for y=1:ny
            V(x,y) = V_temp(y + (x - 1)*ny);
        end
    end
end

% Total Current I 
function I = TotalCur(nx, ny, V1, SigmaConduct, SigmaInsulate, W, L)

    CM = SigmaConduct*ones(nx, ny);
    CM(1:W,(1:L)+ny/2-L/2) = SigmaInsulate;
    CM((1:W)+nx-W,(1:L)+ny/2-L/2) = SigmaInsulate;
    V = NuS(nx, ny, CM, Inf, Inf, 0, V1);
    [Ex, Ey] = gradient(V);
    Ex = -Ex;
    Ey = -Ey;
    Jx = CM.*Ex;
    I = (abs(sum(Jx(1,:))) + abs(sum(Jx(nx,:))))/2;
end