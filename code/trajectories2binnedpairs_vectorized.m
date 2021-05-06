% Dhruv Balwada
% 5 May 2021

% Vectorized calculation of pairs
% Follows the code that was tested in test_new_SFcalc.m

clear all
close all

%%

%traj = load('~/work_root/GoMexico_drifters/GLAD_15min_filtered/traj_mat_GLAD_15min_04_May_2021.mat');

traj = load('~/work_root/GoMexico_drifters/LASER_SPOT_15min_filtered/traj_mat_LASER_15min_04_May_2021.mat');
%% Break up by time and separation bins
% took 58.9s for GLAD
% took 146s for LASER 
tic 
for i=1:length(traj.T_axis)
    id = find(~isnan(traj.trajmat_X(i,:)));
    if mod(i,300)==0
        disp(i)
    end
    
    if length(id)>1
        X = traj.trajmat_X(i,id)';
        Y = traj.trajmat_Y(i,id)';
        U = traj.trajmat_U(i,id)';
        V = traj.trajmat_V(i,id)';
        
        Xvec = [X, Y];
        
        pairs_time(i).dist = pdist(Xvec, @dist_geo);
        
        rx = pdist(Xvec, @dist_rx);
        ry = pdist(Xvec, @dist_ry);
        
        magr = sqrt(rx.^2 + ry.^2);
        
        rx = rx./magr; ry = ry./magr;
        
        dux = pdist(U, @dist_du);
        duy = pdist(V, @dist_du);
        
        pairs_time(i).dul = dux.*rx + duy.*ry;
        pairs_time(i).dut = duy.*rx - dux.*ry;
    else
        pairs_time(i).dul = NaN;
        pairs_time(i).dut = NaN;
        pairs_time(i).dist = NaN;
    end
end
toc

%%
save structure_pairs_LASER.mat pairs_time -v7.3


%%
function rx = dist_rx(XI, XJ)
rx = (XI(:,1) - XJ(:,1)).*cosd(0.5*(XI(:,2) + XJ(:,2)));
end
function ry = dist_ry(XI, XJ)
ry = (XI(:,2) - XJ(:,2));
end
function du = dist_du(UI, UJ)
du = UI - UJ;
end
function dist = dist_geo(XI,XJ)
X = abs(XI(:,1) - XJ(:,1)) .*cosd(0.5*(XI(:,2)+XJ(:,2))) *111321;
Y = abs(XI(:,2) - XJ(:,2)) *111321;
dist = sqrt(X.^2 + Y.^2);
end
