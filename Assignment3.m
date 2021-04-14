%Khaled AbouShaban 101042658
%Assignment 3 for Elec4700, Monte-Carlo/Finite Difference Method

%Part 1 Electron Modelling
%This part is from Assignment 1 but has been eddited for assignment 3
% Monte-Carlo Modeling of Electron Transport without the bottle-neck
%
%% PLEASE LOOK AT THIS PART FOR FIXING MY ASSIGNMENT 1 grades

clear all
close all
set(0,'DefaultFigureWindowStyle','docked')
% Constants

m0 = 9.10938356e-31; %electron mass
k = 1.38064852e-23; %Boltzmann's constant
T = 300; %Temp in Kelvin
m = 0.26*m0; %Effective mass of electron
e = 1.60217662e-19; %electron charge
electron_concentration = 10e15; %electron concentration of 10^15
Volt_applied = 0.1; %
% Nominal size is 200 nm x 100 nm
x_limit = 200e-9;
y_limit = 100e-9;
iters =200;
v_thermal = sqrt(2*k*T/m);
time_step = y_limit/v_thermal/100;
rezo_y = 18;
rezo_x = rezo_y * 2;

% 1 a) Electric Field of Electron b) Force c) Acceleration
voltage = 0.1;
E_field = voltage/x_limit %a
E_force = e*E_field %b
E_Acceleration = E_force/m %c

timestep =y_limit/v_thermal/100;
num_particles = 10000; 
Iter_sim = 1000;
collision_time = 0.2e-12;

% creating arrays for use

particle_array = []; %Stores Particle Array 
temp_array = []; %Stores Temp Array
den_array = []; %Stores Density Array
particle_xpos = []; %Stores X 
particle_ypos = []; %Stores Y
currentx_array = []; %Stores X Current Array
temp_map = []; %Tempreture Map Array
scat_num = 0; 
scat_time = 0;
scat_path = 0;

scat_probability = 1 - exp(-timestep/collision_time);

%initializing the particles

for i = 1:num_particles
    %init positions
    particle_array(i,1) = rand * x_limit;
    particle_array(i,2) = rand * y_limit;
    vth = randn*v_thermal + v_thermal;
    
    %init velocity
    
    particle_array(i,3) = vth; %thermal velocity magnitude
    particle_array(i,4) = ((rand * 2) - 1)*vth; 
    particle_array(i,5) = sqrt(particle_array(i,3)^2 - particle_array(i,4)^2); 
    if(rand > 0.5)
        particle_array(i,5) = particle_array(i,5) * -1;
    end
    
    particle_array(i,6) = 0; %set time=0
        
end

% Starting Simulation
for simCount = 1:Iter_sim
    currTime = simCount * timestep;
        
    % Scattering
    for i = 1:num_particles
        particle_array(i,6) = particle_array(i,6) + timestep;
        
        if(rand <= scat_probability) 
            scat_path = scat_path + (particle_array(i,6)*particle_array(i,3));
            scat_time = scat_time + particle_array(i,6);
            particle_array(i,6) = 0;            
            scat_num = scat_num + 1;
            particle_array(i,3) = randn*v_thermal + v_thermal;
            particle_array(i,4) = ((rand * 2) - 1)*particle_array(i,3); 
            particle_array(i,5) = sqrt(particle_array(i,3)^2 - particle_array(i,4)^2); 
            if(rand > 0.5)
                particle_array(i,5) = particle_array(i,5) * -1;
            end
        end
    end
    
   
    particle_array(:,4) = particle_array(:,4) + timestep * E_Acceleration;
    
    % new particle positions
    new_xPos = particle_array(:,1) + timestep * particle_array(:,4);
    new_yPos = particle_array(:,2) + timestep * particle_array(:,5);    
    
    
    %checking boundary conditions
    for i = 1:num_particles
        if(new_xPos(i) < 0) 
            new_xPos(i) = new_xPos(i)+x_limit;
        elseif(new_xPos(i) > x_limit) 
            new_xPos(i) = new_xPos(i) - x_limit;
        end
        
        if(new_yPos(i) < 0) 
            new_yPos(i) = abs(new_yPos(i));
            particle_array(i,5) = particle_array(i,5) * -1; 
        elseif(new_yPos(i) > y_limit) 
            new_yPos(i) = 2*y_limit - new_yPos(i);
            particle_array(i,5) = particle_array(i,5) * -1; 
        end
       
    end
    
    % new pos
    particle_array(:,1) = new_xPos;
    particle_array(:,2) = new_yPos;
    


    % save pos 
    particle_pos_x_arr(1:7,simCount) = particle_array(1:7,1);
    particle_pos_y_arr(1:7,simCount) = particle_array(1:7,2);

  
    % storing the X current
    xCurrent_Array(simCount) = e * electron_concentration * mean(particle_array(:,4)) * x_limit;


    
    if(simCount == Iter_sim) 
        
        x_d = x_limit / rezo_x;
        y_d = y_limit / rezo_y;
        
        for yCount = 1:rezo_y
            for xCount = 1:rezo_x
                densityMap_Array(xCount,yCount) = 0;
                temp_maparr(xCount,yCount) = 0;
                
                for count_2 = 1:num_particles
                    if(particle_array(count_2,1) >= (xCount*x_d-x_d) && particle_array(count_2,1) < (xCount*x_d))
                        if(particle_array(count_2,2) >= (yCount*y_d-y_d) && particle_array(count_2,2) < (yCount*y_d))
                            densityMap_Array(xCount,yCount) = densityMap_Array(xCount,yCount) + 1;
                            temptemp = sqrt(particle_array(count_2,4).^2+particle_array(count_2,5).^2)*m/(3*k);
                            temp_maparr(xCount,yCount) = temp_maparr(xCount,yCount) + temptemp;
                        end                        
                    end
                end                
            end
        end   
        temp_maparr = temp_maparr ./ densityMap_Array;
    end
end


%plot of the trajectory of the particles
figure(1);
hold on;
for i = 1:7
    plot(particle_pos_x_arr(i,:),particle_pos_y_arr(i,:));
end
% xlim([0,x_limit]);
% ylim([0,y_limit]);
title(['Particle Trajectory, Simulation Count ', num2str(simCount), ' Timestep: ', num2str(timestep)]);
xlabel('X (m)');
ylabel('Y (m)');

%Plot Current/Time
figure(2);

x = linspace(timestep,Iter_sim*timestep,Iter_sim);

plot(x,xCurrent_Array);
title('X Current over Time');
xlabel('Time (s)');
ylabel('Current (A)');

%plot the density map
figure(3);
[x_mesh, y_mesh] = meshgrid(linspace(0,x_limit,rezo_x),linspace(0,y_limit,rezo_y));
surf(x_mesh,y_mesh,transpose(densityMap_Array));
title('Density Map');
xlabel('x (m)');
ylabel('y (m)');
grid on;
colorbar ;

% figure(4);
% surf(x_mesh,y_mesh);
% colorbar;
% grid on;
% title('Density Map')
% xlabel('x position (nm)');
% ylabel('y position (nm)');
% colorbar;
% hold off;

%plotting the temp map


figure(4);
[x_mesh, y_mesh] = meshgrid(linspace(0,x_limit,rezo_x),linspace(0,y_limit,rezo_y));
surf(x_mesh,y_mesh,transpose(temp_maparr));
title('Temperature Map');
xlabel('x (m)');
ylabel('y (m)');
grid on;
colorbar ;

%% 

% 
% Part 2 using FDM from assignment 2 to calc the elctric field inluding the Bottle-Neck Sim
% This part uses Finite Difference Method in Assignment-2 to calculate the electric field
% and then provide a field for the Monte-Carlo bottle-neck simulation.
% 
% PLEASE LOOK AT THIS PART FOR FIXING MY ASSIGNMENT 2 grades


ySize = 100; %size of the matrix
xSize = ySize*2; %size of the matrix
xDel = x_limit/xSize; %delta X
yDel = y_limit/ySize; %delta Y
maxIterations = 10000; %maximum number of iterations

v_matrix = zeros(ySize,xSize); %create our matrix of values
oldv_matrix = zeros(ySize,xSize); %another matrix to be used to clone the new one
vecMatrix = zeros(ySize,xSize); %matrix of the vectors
[x_mesh, y_mesh] = meshgrid(linspace(0,x_limit,xSize),linspace(0,y_limit,ySize));

% Setting Boundary 
left = 0.1;
right = 0;
north = 0; 
south = 0; 

y_lim2 = 25;
x_lim2 = 35;

for i = 1:maxIterations %run till max reached
    for m = 1:xSize % run through x 
        for n = 1:ySize %run through y
            
            if(m == 1) %left bound
                v_matrix(n,m) = left;
            elseif(m == xSize) %right bound
                v_matrix(n,m) = right;
            elseif(n == 1) %North bound
                v_matrix(n,m) = (oldv_matrix(n,m-1) + oldv_matrix(n,m+1) + oldv_matrix(n+1,m))/3;
            elseif(n == ySize) %South bound
                v_matrix(n,m) = (oldv_matrix(n,m-1) + oldv_matrix(n,m+1) + oldv_matrix(n-1,m))/3;
            else
                v_matrix(n,m) = (oldv_matrix(n,m-1) + oldv_matrix(n,m+1) + oldv_matrix(n-1,m) + oldv_matrix(n+1,m))/4;
            end   
            
            if(m>=xSize/2-x_lim2/2 && m<=xSize/2+x_lim2/2 && (n<=y_lim2 || n>=ySize-y_lim2)) %inside one of the boxes
                v_matrix(n,m) = v_matrix(n,m)*10^-2;
            end
        end
    end        
    [xVectors, yVectors] = gradient(oldv_matrix(:,:));   
  
    oldv_matrix = v_matrix;
end

figure(5)
%CAP = x_mesh/y_mesh,
surf(x_mesh, y_mesh, v_matrix); %surf ( v_matrix ) ;
shading interp
title('V(x,y) for Bottle-Neck');
xlabel('x (m)');
ylabel('y (m)');
grid on;  
hold on;

figure(6)
quiver(x_mesh, y_mesh);
xlim([0 100]);
ylim([0 100]);
title('Electric Field for Bottle-Neck');
xlabel('x (m)');
ylabel('y (m)');
grid on;

particle_array = [];

tempArray = []; %array storing temp
particle_pos_x_arr = []; %array for x array at the 7th pos
particle_pos_y_arr = []; %aarray for y array at the 7th pos
scat_path = 0; 
scat_time = 0; 
scat_num = 0; 
xCurrent_Array = []; 
electron_concentration = 10e15; 
densityMap_Array = [];
temp_maparr = [];
rezo_y = 20;
rezo_x = rezo_y * 2;

scat_probability = 1 - exp(-timestep/collision_time); %probability of particle getting scattered


%initiazation of particles
num_particles = 1000;
for i = 1:num_particles
    
    particle_array(i,1) = rand * x_limit;
    particle_array(i,2) = rand * y_limit;
    vth = randn*v_thermal + v_thermal;
    
    while(particle_array(i,1)>=(x_limit/2-x_lim2/2*1e-9) && particle_array(i,1)<=(x_limit/2+x_lim2/2*1e-9) && (particle_array(i,2)<=(y_lim2*1e-9) || particle_array(i,2)>=(y_limit-y_lim2*1e-9)))
        %choose other location
        particle_array(i,1) = rand * x_limit;
        particle_array(i,2) = rand * y_limit;
    end
    
    %initate velocity
    particle_array(i,3) = vth; % setting vth
    particle_array(i,4) = ((rand * 2) - 1)*vth; 
    particle_array(i,5) = sqrt(particle_array(i,3)^2 - particle_array(i,4)^2); 
    if(rand > 0.5)
        particle_array(i,5) = particle_array(i,5) * -1;
    end
    
    particle_array(i,6) = 0; %set time since last scatter to 0
        
end

%run simulation
Iter_sim = 1000;
for simCount = 1:Iter_sim
    currTime = simCount * timestep;
        
    %scatter particles
    for i = 1:num_particles
        %update time since last scatter
        particle_array(i,6) = particle_array(i,6) + timestep;
        
        if(rand <= scat_probability) %scatter
            scat_path = scat_path + (particle_array(i,6)*particle_array(i,3)); 
            scat_time = scat_time + particle_array(i,6); 
            particle_array(i,6) = 0; %reset time          
            scat_num = scat_num + 1;
            particle_array(i,3) = randn*v_thermal + v_thermal; 
            particle_array(i,4) = ((rand * 2) - 1)*particle_array(i,3); 
            particle_array(i,5) = sqrt(particle_array(i,3)^2 - particle_array(i,4)^2); 
            if(rand > 0.5)
                particle_array(i,5) = particle_array(i,5) * -1;
            end
        end
    end
      
    for count_2 = 1:num_particles
        for yCount = 1:ySize
            for xCount = 1:xSize            
                if(particle_array(count_2,1) >= (xCount*xDel-xDel) && particle_array(count_2,1) < (xCount*xDel))
                    if(particle_array(count_2,2) >= (yCount*yDel-yDel) && particle_array(count_2,2) < (yCount*yDel))
                        particle_array(count_2,4) = particle_array(count_2,4) + timestep * -xVectors(yCount,xCount)*e/m/xDel;
                        particle_array(count_2,5) = particle_array(count_2,5) + timestep * -yVectors(yCount,xCount)*e/m/xDel; 
                    end                        
                end
            end                
        end
    end  
    
    %update pos of particles
    new_xPos = particle_array(:,1) + timestep * particle_array(:,4);
    new_yPos = particle_array(:,2) + timestep * particle_array(:,5);    
    
    
    %checking boundary conditions
    for i = 1:num_particles
        if(new_xPos(i) < 0) 
            new_xPos(i) = new_xPos(i)+x_limit;
            
        elseif(new_xPos(i) > x_limit) 
            new_xPos(i) = new_xPos(i) - x_limit;
        end
        
        if(new_yPos(i) < 0) 
            new_yPos(i) = abs(new_yPos(i));
            particle_array(i,5) = particle_array(i,5) * -1; 
        elseif(new_yPos(i) > y_limit) 
            new_yPos(i) = 2*y_limit - new_yPos(i);
            particle_array(i,5) = particle_array(i,5) * -1; 
        end
       
        % if in the box
        bx_leftSide = x_limit/2-x_lim2/2*1e-9;
        bx_rightSide = x_limit/2+x_lim2/2*1e-9;
        box_bottom = y_lim2*1e-9;
        box_top = y_limit-y_lim2*1e-9;
        
        if(new_xPos(i) >= bx_leftSide && new_xPos(i) <= bx_rightSide && (new_yPos(i) <= box_bottom || new_yPos(i) >= box_top))
            if(particle_array(i,1) < bx_leftSide)
                tempx = bx_leftSide - abs(new_xPos(i) - particle_array(i,1));
                if(~(tempx >= bx_leftSide && tempx <= bx_rightSide && (new_yPos(i) <= box_bottom || new_yPos(i) >= box_top))) %not in box
                    new_xPos(i) = tempx;
                    particle_array(i,4) = particle_array(i,4) * -1;
                end
            elseif(particle_array(i,1) > bx_rightSide)
                tempx = bx_rightSide + abs(new_xPos(i) - particle_array(i,1));
                if(~(tempx >= bx_leftSide && tempx <= bx_rightSide && (new_yPos(i) <= box_bottom || new_yPos(i) >= box_top))) %not in box
                    new_xPos(i) = tempx;
                    particle_array(i,4) = particle_array(i,4) * -1; 
                end
            elseif(particle_array(i,2) < box_top)
                
                tempy = box_top - abs(new_yPos(i) - particle_array(i,2));
                if(~(new_xPos(i) >= bx_leftSide && new_xPos(i) <= bx_rightSide && (tempy <= box_bottom || tempy >= box_top))) %not in box
                    new_yPos(i) = tempy;
                    
                    particle_array(i,5) = particle_array(i,5) * -1;
                end
            elseif(particle_array(i,2) > box_bottom)
                tempy = box_bottom + abs(new_yPos(i) - particle_array(i,2));
                if(~(new_xPos(i) >= bx_leftSide && new_xPos(i) <= bx_rightSide && (tempy <= box_bottom || tempy >= box_top))) %not in box
                    new_yPos(i) = tempy;
                    particle_array(i,5) = particle_array(i,5) * -1; 
                end
            end
        end
        
    end
    
   
    
    %setting new positions
    particle_array(:,1) = new_xPos;
    particle_array(:,2) = new_yPos;
    
  
    %save the particle position for plotting
    particle_pos_x_arr(1:7,simCount) = particle_array(1:7,1);
    particle_pos_y_arr(1:7,simCount) = particle_array(1:7,2);

  

    %store the X current
    xCurrent_Array(simCount) = e * electron_concentration * mean(particle_array(:,4)) * x_limit;

    pause(0.0001); 
    
    if(simCount == Iter_sim) %last iter
        
        x_d = x_limit / rezo_x;
        y_d = y_limit / rezo_y;
        
        for yCount = 1:rezo_y
            for xCount = 1:rezo_x
                densityMap_Array(xCount,yCount) = 0;
                temp_maparr(xCount,yCount) = 0;
                
                for count_2 = 1:num_particles
                    if(particle_array(count_2,1) >= (xCount*x_d-x_d) && particle_array(count_2,1) < (xCount*x_d))
                        if(particle_array(count_2,2) >= (yCount*y_d-y_d) && particle_array(count_2,2) < (yCount*y_d))
                            densityMap_Array(xCount,yCount) = densityMap_Array(xCount,yCount) + 1;
                            temptemp = sqrt(particle_array(count_2,4).^2+particle_array(count_2,5).^2)*m/(3*k);
                            temp_maparr(xCount,yCount) = temp_maparr(xCount,yCount) + temptemp;
                        end                        
                    end
                end                
            end
        end   
        temp_maparr = temp_maparr ./ densityMap_Array;
        temp_maparr(isnan(temp_maparr)) = 0; %detemine which element are Nan
    end
end


%plot the trajectory of the particles
figure(7);
hold on;
for i = 1:7
    plot(particle_pos_x_arr(i,:),particle_pos_y_arr(i,:));
end
xlim([0,x_limit]);
ylim([0,y_limit]);
title(['Particle Trajectory, Simulation Count: ', num2str(simCount), ', Timestep: ', num2str(timestep)]);
xlabel('X (m)');
ylabel('Y (m)');

%Plot the current across X
figure(8);
x = linspace(timestep,Iter_sim*timestep,Iter_sim);
plot(x,xCurrent_Array);
title('Current/Time Part 3');
xlabel('time (s)');
ylabel('current (A)');

%Plot the density map
figure(9);
[x_mesh, y_mesh] = meshgrid(linspace(0,x_limit,rezo_x),linspace(0,y_limit,rezo_y));
surf(x_mesh,y_mesh,transpose(densityMap_Array));
title('Density Map Plot');
xlabel('x (m)');
ylabel('y (m)');
grid on;
colorbar;

%plotting the temperature map
figure(10);
[x_mesh, y_mesh] = meshgrid(linspace(0,x_limit,rezo_x),linspace(0,y_limit,rezo_y));
surf(x_mesh,y_mesh,transpose(temp_maparr));
title('Temperature Map Plot');
xlabel('x (m)');
ylabel('y (m)');
grid on;
colorbar

avgCurrent = mean(xCurrent_Array(100:end)); %Average Current
BottleneckWidth = 100-y_lim2*2;  %Bottle-Neck Width


%plot the avg current vs bottleneck width

arr_average = [0.053121, 0.069143, 0.088485, 0.096499, 0.10678, 0.11388, 0.1328, 0.12328, 0.12880];
bott_Width = [10, 20, 30, 40, 50, 60, 70, 80, 90]; %Seeting width arr

figure(11);
plot(bott_Width,arr_average);
title('Average Current Vs Bottle-neck Width');
xlabel('Bottleneck Width (nm)');
ylabel('Current (A)');