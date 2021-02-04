%% Trajectories to binned pairs
% Dhruv Balwada, 29 Jan 2021
%
% The first step is to go from the trajectory data to
% pair velocities, which are binned in separation bins. 
%
% The trajectory data was downloaded from: 
% https://data.gulfresearchinitiative.org/metadata/R1.x134.073:0004
% (it is not in the repo because it is a large file ~ 30mb)

clear all
close all

traj = load ('/Users/dhruvbalwada/Google Drive/Main_projects/GOM/GLAD_data/glad_traj_297.mat');

%% Calculate time series of pair separation from trajectories

sep = calculate_seperation_timeseries(traj);
% this operation is quite slow at this point

%% Define separation bins

gamma = 1.5;

dist_bin(1) = 0.1; % in m
dist_bin = gamma.^[0:100]*dist_bin(1);

% get rid of bins that are too large
id = find(dist_bin>800*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

%% Bin the separation data 

% loop for different distance classes
% this is quite slow.
for i =1:length(dist_axis)
    dull = []; dutt = [];
    % loop over different pairs
    for j = 1:length(sep)
        dull_temp = []; dutt_temp = [];
        % find id of the pairs in a particular geographical regime (within a certain distance from each other)
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i));
        
        % loop over the different pairs that lie in the range
        for k =1:length(id)
            % components of vector joining the two particles
            rx(k) = (sep(j).X1(id(k)) - sep(j).X2(id(k)))*cosd(0.5*(sep(j).Y1(id(k))+sep(j).Y2(id(k))));
            ry(k) = (sep(j).Y1(id(k)) - sep(j).Y2(id(k)));
            magr(k) = sqrt(rx(k).^2+ry(k).^2);
            
            % normalize to unit vectors
            rx(k) = rx(k)/magr(k); ry(k) = ry(k)/magr(k);
            
            % components of velocity differences
            dux(k) = (sep(j).U1(id(k))-sep(j).U2(id(k)));
            duy(k) = (sep(j).V1(id(k))-sep(j).V2(id(k)));
            
            % convert to longitudnal and transverse structure functions
            dull_temp(k) = dux(k)*rx(k) + duy(k)*ry(k);
            dutt_temp(k) = duy(k)*rx(k) - dux(k)*ry(k);

        end
        if ~isempty(dull_temp)
            dull = [dull; dull_temp'];
            dutt = [dutt; dutt_temp'];
        end
    end
    
    struct_pairs(i).dull = dull;
    struct_pairs(i).dutt = dutt;
    
    disp(i)
end

%% Save the pair in different separation bins for use later

save structure_pairs.mat struct_pairs dist_axis dist_bin -v7.3
% this is a large file ~2.4gb (this is also not in the repo)