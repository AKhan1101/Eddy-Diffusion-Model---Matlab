function conc = point_eddy_diff_conc(pos,k,t,lambda,mass,room,source)
% 16-11-23 - Ahmed Khan
% This function calculates the concentration of a pollutant at a
% discretised point in space and time using the eddy diffusion model. 
%% Inputs
% node_pos - a 3x1 [x,y,z] vector of nodal co-ordinates.
% k - eddy diffusion coefficient.
% t - simualtion time [s]
% lambda - fresh air change rate [kg/min] 
% mass - 
% room - [length, width, height] of the room.
% source - [x,y,z] cartesian co-ordinates. 
% img_source - [x_1, y_1, z_1; x_2, y_2, z_2;...] Average Co ordinated of
% each wall 
%%
ACH = lambda; 
q = mass; 
wall_coeff = zeros(6,1);

w = room(1);
l = room(2); 
h = room(3); 

%% Main loop

for n = 1:6                                                                % Need to make so that each face/img source has the correct co-ordinated called. 
    a = 0:6;

    rx(n) = exp((-(pos(1)+2*a(n)*w-source(1))^2)./4*k*t)+...
        exp((-(pos(1)+2*a(n)*w+source(1))^2)./4*k*t);
    ry(n) = exp((-(pos(2)+2*a(n)*l-source(2))^2)./4*k*t)+...
        exp((-(pos(2)+2*a(n)*l+source(2))^2)./4*k*t); 
    rz(n) = exp((-(pos(3)+2*a(n)*h-source(3))^2)./4*k*t)+...
        exp((-(pos(3)+2*a(n)*h+source(3))^2)./4*k*t); 

    wall_coeff(n) = rx(n).*ry(n).*rz(n);

end

wall_coeff_sum = sum(wall_coeff);

conc = ((mass*exp(-lambda.*t))./(8*(pi*k*t)^3/2))*  wall_coeff_sum;   
end
