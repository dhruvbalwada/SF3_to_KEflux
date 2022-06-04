% Dhruv Balwada
% 4 May 2021

% Function to recast drifter structure format into arrays with shape 
% (num_time_steps, num_trajectories). 

clear all
close all

% open structure file
sel_data = 'LASER'; 

if sel_data == 'GLAD'
   load ../GLAD_15min_filtered/traj_struct_GLAD_15min_03_May_2021.mat
elseif sel_data == 'LASER'
    load ../LASER_SPOT_15min_filtered/traj_structs_LASER_15min_03_May_2021.mat
    drifter=drifterL; % select drogued for LASER
end


%%
ndrifts = length(drifter);

tmin = zeros(ndrifts,1);
tmax = zeros(ndrifts,1);
dt_mean = zeros(ndrifts,1);
dt_max = zeros(ndrifts,1);
dt_min = zeros(ndrifts,1);

for i = 1:ndrifts
    tmin(i) = min(drifter(i).time); 
    tmax(i) = max(drifter(i).time); 
    
    % check if there are any data gaps (done visually)
    dt_mean(i) = mean(diff(drifter(i).time));
    dt_min(i) = min(diff(drifter(i).time));
    dt_max(i) = max(diff(drifter(i).time));
end

tmin_min = min(tmin);
tmax_max = max(tmax);


%% 
figure,
scatter(dt_min, dt_max) % should be close to 0.0104 
%% 

T_axis = [tmin_min:mean(dt_mean):tmax_max]; 
% This is only approximate, and likely works only upto the nearest minute. 
% However, since the data does not have much variance at scales smaller
% than an hour, it is likely not important. 

%% 
trajmat_X = NaN*ones(length(T_axis), ndrifts);
trajmat_Y = NaN*ones(length(T_axis), ndrifts);
trajmat_U = NaN*ones(length(T_axis), ndrifts);
trajmat_V = NaN*ones(length(T_axis), ndrifts);

for i = 1:ndrifts
    tinit = drifter(i).time(1); 
    tfinal = drifter(i).time(end); 
    
    id_start = find(abs(T_axis-tinit)<6e-4);
    id_end = find(abs(T_axis-tfinal)<6e-4);
    
    len = id_end - id_start + 1; 
    
    trajmat_X(id_start:id_end,i) = drifter(i).lon;
    trajmat_Y(id_start:id_end,i) = drifter(i).lat;
    trajmat_U(id_start:id_end,i) = drifter(i).u;
    trajmat_V(id_start:id_end,i) = drifter(i).v;
end

%%

if sel_data == 'GLAD'
    save ../GLAD_15min_filtered/traj_mat_GLAD_15min_04_May_2021.mat T_axis trajmat_X trajmat_Y trajmat_U trajmat_V
elseif sel_data == 'LASER'
    save ../LASER_SPOT_15min_filtered/traj_mat_LASER_15min_04_May_2021.mat T_axis trajmat_X trajmat_Y trajmat_U trajmat_V
end


%% some rough plots

%
figure
subplot(211)
plot(trajmat_X)

subplot(212)
plot(trajmat_Y)

%% 
figure
subplot(211)
plot(trajmat_U)

subplot(212)
plot(trajmat_V)

%%
figure
plot(trajmat_X, trajmat_Y)
