%% Data distribution
% Dhruv Balwada, 23 December 2021

% Plots of distribution as:
% - Time
% - Separation
% - Space (binning)
%
clear all
close all
%% Load the file where the binned pairs are
GLAD_pairs = load('../data/structure_pairs_GLAD_deep_500m.mat');
LASER_pairs = load('../data/structure_pairs_LASER_deep_500m_box_constrained.mat');

% Load num pairs per bin 
GLAD_numpairs = load('../data/GLAD_npairs_deep500_block_boot_strap_Ldof.mat');
LASER_numpairs = load('../data/LASER_npairs_deep500_box_constrained_block_boot_strap_Ldof.mat');

%% Load trajectory data
GLAD_traj = load('/Users/dhruvbalwada/OneDrive/GoMexico_drifters_data/GLAD_15min_filtered/traj_mat_GLAD_15min_04_May_2021.mat');
temp =  load('../data/traj_depth_GLAD.mat');
GLAD_traj.trajmat_H = temp.Htraj; 

LASER_traj = load('/Users/dhruvbalwada/OneDrive/GoMexico_drifters_data/LASER_SPOT_15min_filtered/traj_mat_LASER_15min_04_May_2021.mat');
temp =  load('../data/traj_depth_LASER.mat');
LASER_traj.trajmat_H = temp.Htraj; 

%% Time plot

GLAD_pairs.npairs  = count_pairs(GLAD_pairs);
LASER_pairs.npairs = count_pairs(LASER_pairs);

%%
cols = colororder;

%% Plots vs time

figure
plot(GLAD_traj.T_axis, GLAD_pairs.npairs, 'linewidth',2, 'color', cols(7,:))
datetick('x','dd-mmm,yy')
set(gca, 'YScale', 'log')
ylim([500 3e5])

set(gca,'fontsize',16, 'fontname','Times')
ylabel('Number of pairs')
xlabel('Date')
print(['./fig_data_dist/dist_time_GLAD.png'],'-dpng', '-r400')

%
figure
plot(LASER_traj.T_axis, LASER_pairs.npairs, 'linewidth',2, 'color', cols(1,:))
datetick('x','dd-mmm,yy')
set(gca, 'YScale', 'log')
ylim([500 3e5])

set(gca,'fontsize',16, 'fontname','Times')
ylabel('Number of pairs')
xlabel('Date')
print(['./fig_data_dist/dist_time_LASER.png'],'-dpng', '-r400')

%% Pairs per sep bin

figure 
loglog(GLAD_numpairs.dist_axis/1e3, GLAD_numpairs.npairs_sep, 'linewidth',2, 'color', cols(7,:))
hold all
loglog(LASER_numpairs.dist_axis/1e3, LASER_numpairs.npairs_sep, 'linewidth',2, 'color', cols(1,:))

legend('GLAD', 'LASER', 'location','northwest')
axis([1e-2 1e3 5e3 1e8])
set(gca,'fontsize',16, 'fontname','Times')
ylabel('Number of pairs')
xlabel('r [km]')

print(['./fig_data_dist/dist_sep.png'],'-dpng', '-r400')

%%

ybin = 24:0.25:31; 
xbin = -91:0.25:-84; 
yaxis = 0.5*(ybin(1:end-1) + ybin(2:end));
xaxis = 0.5*(xbin(1:end-1) + xbin(2:end));

GLAD_nperbin = zeros(length(yaxis), length(xaxis)); 
LASER_nperbin = zeros(length(yaxis), length(xaxis)); 


for i = 1:length(yaxis)
    for j = 1:length(xaxis)
        GLAD_nperbin(i, j) = length(find(GLAD_traj.trajmat_X>= xbin(j) & GLAD_traj.trajmat_X<xbin(j+1) ...
                            & GLAD_traj.trajmat_Y>=ybin(i) & GLAD_traj.trajmat_Y<ybin(i+1)  ...
                            &  GLAD_traj.trajmat_H<-500)); 
        LASER_nperbin(i, j) = length(find(LASER_traj.trajmat_X>= xbin(j) & LASER_traj.trajmat_X<xbin(j+1) ...
                            & LASER_traj.trajmat_Y>=ybin(i) & LASER_traj.trajmat_Y<ybin(i+1) & ...
                            LASER_traj.trajmat_H<-500));                         
    end
end


%% Full domain figure 

figure
m_proj('Lambert','lon',[-93 -82],'lat',[20 32]); 


m_pcolor(xaxis, yaxis, log10(GLAD_nperbin))
colormap('hot');
brighten(.5);
colorbar

caxis([0.5 5])
hold all

m_grid('linestyle','none','tickdir','in','linewidth',1,'FontSize',20,'FontName','Times','interpreter','latex');


m_coast('patch',[0.5 0.5 0.5])
set(gca,'fontsize',16, 'fontname','Times')

title('GLAD')
print(['./fig_data_dist/bins_glad.png'],'-dpng', '-r400')

%%
figure
m_proj('Lambert','lon',[-93 -82],'lat',[20 32]); 


m_pcolor(xaxis, yaxis, log10(LASER_nperbin))
colormap('hot');
brighten(.5);
colorbar

caxis([0.5 5])
hold all

m_grid('linestyle','none','tickdir','in','linewidth',1,'FontSize',20,'FontName','Times','interpreter','latex');


m_coast('patch',[0.5 0.5 0.5])
set(gca,'fontsize',16, 'fontname','Times')

title('LASER')
print(['./fig_data_dist/bins_laser.png'],'-dpng', '-r400')

%%
function [npairs] = count_pairs(ds)
    npairs = zeros(length(ds.pairs_time),1);
    for i = 1:length(ds.pairs_time)    
        npairs(i) = length(find(~isnan(ds.pairs_time(i).dul)));
    end
end