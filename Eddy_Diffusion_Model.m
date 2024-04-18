%% Inputs 
room_width = 10;                                                           % [m]
room_length = 10;                                                          % [m] 
room_height = 5;                                                           % [m]

% room_length = 27.6;
% room_width = 13.7;
% room_height = 4.5;

nNodes_x = 125;                                                          % Number of points within the room. 
nNodes_y = 125;                                                          % Number of points within the room. 
nNodes_z = 9;                                                            % Number of points within the room.

total_time = 10;                                                           % time - [seconds]
time_step = 2;                                                             

k = 2E-03;

mass = 1E-03;

lambda = 0.005;

source = [0,0,0];
%% Pre Processing

nodes_x = linspace (-room_width./2,room_width./2,nNodes_x);              % Discretising node co-ordinates
nodes_y = linspace (-room_length./2,room_length./2,nNodes_y);            %
nodes_z = linspace (0,room_height,nNodes_z);                             %
            
ntime = (total_time./time_step);                                           % number of time steps.

conc = zeros(nNodes_x,nNodes_y,nNodes_z,ntime);                            % Initialising grid
                                                                        
% delta t - [seconds] 

source_edy = [source(2),source(1),source (3)];                             % Cartesian co-ordinates [m] 

room = [room_width, room_length, room_height];                             %

t_point = linspace(time_step,total_time,ntime);                            % Time points solution will be calculated @ 
                         
%% Calculating the concentration
% Formulation from:

%Input of time is incorrect. 
% Release co-ordionates may be inputted incorrectly. 
delta_conc = zeros (nNodes_x,nNodes_y,nNodes_z,length(t_point));

for t = 1:length(t_point)                                                  % Loop trough time 
    for x = 1:nNodes_x                                                     % Loop trough x co-ordinates 
        for y = 1:nNodes_y                                                 % Loop trough y co-ordinates
            for z = 1:nNodes_z                                             % Loop trough z co-ordinates
                node_pos = [nodes_x(x); nodes_y(y); nodes_z(z)];           % Vectorising the position solving for
                delta_conc(x,y,z,t) = point_eddy_diff_conc_validation(node_pos,k,...  % Calculating the point concentration. 
                    t_point(t),lambda,mass,room,source_edy);
 
            end
        end
    end
end

for x = 1:nNodes_x                                                         % Integrating
    for y = 1:nNodes_y
        for z = 1:nNodes_z

            for t = 1:length(t_point)

                delta_conc_sum = reshape(delta_conc(x,y,z,:),1,length(t_point));                
                conc (x,y,z,:) = cumtrapz(delta_conc_sum,t_point);

            end
            conc_Tavg (x,y,z) = mean (delta_conc_sum); 

        end
    end
end

%% Plotting
plot = true;

if plot == true 
    for i = 1:ntime                                                            % Plotting. 
    
        figure ()
        hold on
        [a,b] = contourf(delta_conc(:,:,1,i));
        % colourmap(delta_conc(:,:,1,i))
        xticks = nodes_x;
        yticks = nodes_y; 
        colorbar
        title ('T = ', t_point(i))
        
    end
    hold off 
end