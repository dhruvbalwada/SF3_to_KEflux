% Dhruv Balwada
% 16 May 2021

% Estimate the depth at the location of the drifters 

clear all
close all

%%

traj = load('~/work_root/GoMexico_drifters/GLAD_15min_filtered/traj_mat_GLAD_15min_04_May_2021.mat');

%traj = load('~/work_root/GoMexico_drifters/LASER_SPOT_15min_filtered/traj_mat_LASER_15min_04_May_2021.mat');

%% Add check of depth 
% We want to select only the drifters in deep waters

% Get bathy
[elev, latg, long] = mygrid_sand2([-99 -79 17 33]);
long = long-360; 

%% Calculate the depth at the location of the trajectory 

Htraj = NaN*ones(size(traj.trajmat_X));
%
for i = 1:size(traj.trajmat_X,1) 
    for j = 1:size(traj.trajmat_X,2)
        X = traj.trajmat_X(i, j);
        Y = traj.trajmat_Y(i, j);
        
        if ~isnan(X)
            idx = find(long<= X, 1, 'last');
            idy = find(latg<= Y, 1, 'last'); 

            Htraj(i,j) = elev(idy, idx); 
        end
    end
end


%%
save ../data/traj_depth_GLAD.mat Htraj
%save ../data/traj_depth_LASER.mat Htraj

%%
figure
subplot(211)
contourf(long, latg, elev)
hold all
plot(traj.trajmat_X(:,1), traj.trajmat_Y(:,1))

subplot(212)
plot(Htraj(:,1))
%%
figure 
hist(Htraj(:))