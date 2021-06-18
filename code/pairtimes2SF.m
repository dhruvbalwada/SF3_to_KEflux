% Dhruv Balwada
% 6 May 2021

% Here we estimate structure functions from the pairs that are formed as a
% function of time (using trajectories2binnedpairs_vectorized.m).

clear all
close all

%load('../data/structure_pairs_GLAD.mat');
%load('../data/structure_pairs_LASER.mat');
load('../data/structure_pairs_GLAD_deep_500m.mat');
%load('../data/structure_pairs_LASER_deep_500m.mat');
%load('../data/structure_pairs_LASER_deep_500m_box_constrained.mat');

%%
traj = load('/Users/dhruvbalwada/work_root/GoMexico_drifters/GLAD_15min_filtered/traj_mat_GLAD_15min_04_May_2021.mat');
%traj = load('/Users/dhruvbalwada/work_root/GoMexico_drifters/LASER_SPOT_15min_filtered/traj_mat_LASER_15min_04_May_2021.mat');
%%

tpts = length(pairs_time);
%
npairs = zeros(tpts,1);
for i = 1:tpts
    npairs(i) = length(find(~isnan(pairs_time(i).dul)));
end

%%
close all 
figure
plot(traj.T_axis, npairs, 'linewidth',2)

%T_ticks = datenum(['21-Jul-2012'; '21-Aug-2012'; '21-Sep-2012'; '21-Oct-2012']);
T_ticks = datenum(['21-Jan-2016'; '21-Feb-2016'; '21-Mar-2016'; '21-Apr-2016']);

xtick(T_ticks)
%xticklabels(['21-Jul-2012'; '21-Aug-2012'; '21-Sep-2012'; '21-Oct-2012'])
xticklabels(['21-Jan-2016'; '21-Feb-2016'; '21-Mar-2016'; '21-Apr-2016'])
title('Num Pairs LASER')
%title('Num Pairs GLAD')
set(gca,'fontsize',16)
%print('./figures/GLAD_num_pairs.eps','-depsc', '-r400')
%print('./figures/LASER_num_pairs.eps','-depsc', '-r400')
print('./figures/LASER_num_pairs_box.eps','-depsc', '-r400')

%% % Align pairs in a single vector
% 
dul = zeros(sum(npairs),1);
dut = zeros(sum(npairs),1);
dist = zeros(sum(npairs),1);

%
% estimate num of pairs

empty1 = 1;
for i = 1:tpts % time loop
    if npairs(i) == 0
        continue
    end
    
    if npairs(i) == 1
        dist(empty1) = pairs_time(i).dist;
        dul(empty1) = pairs_time(i).dul;
        dut(empty1) = pairs_time(i).dut;
        empty1 = empty1+1;
    end
    
    if npairs(i) >1
        dist(empty1: empty1+npairs(i)-1) = pairs_time(i).dist;
        dul(empty1: empty1+npairs(i)-1) = pairs_time(i).dul;
        dut(empty1: empty1+npairs(i)-1) = pairs_time(i).dut;
        empty1 = empty1+npairs(i);
    end
    
end

%%
clear pairs_time

%%  Compute Statistics in vel-diff vs separation space

% generate distance axis

gamma = 1.5;

dist_bin(1) = 10; % in m
dist_bin = gamma.^[0:100]*dist_bin(1);
id = find(dist_bin>1000*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));


% Generate vel axis
vel_bins = linspace(-2, 2, 50);
vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));
%% Compute the histogram
edges = {dist_bin vel_bins};
%
Ndul = hist3([dist, dul], 'Edges', edges);

% Normalize histogram per dist_bins
Ndul_norm = Ndul./sum(Ndul,2);

%%
% This figure looks nice. The one downside is that that the variance
% increases with dist_bins, and so there is not much resolution at the
% smallest scales. 
% One way to make this plot better would be to normalize by variance. 
figure,
contourf(log10(dist_axis/1e3), vel_axis, log(Ndul_norm(1:end-1, 1:end-1))')
grid on
c = colorbar;
c.Label.String = 'log_{10}(PDF)';
xticks(log10([0.01 0.1 1 10 100 1000]))
xticklabels([0.01 0.1 1 10 100 1000])
xlabel('r [km]')
ylabel('\delta u_l [m/s]')
%set(h,'fontsize',16)
set(gca,'fontsize',16, 'fontname','Times')

%title('GLAD')
%print('./figures/GLAD_PDF.eps','-depsc', '-r400')
%print('./figures/LASER_PDF.eps','-depsc', '-r400')
print('./figures/LASER_PDF_box.eps','-depsc', '-r400')
%% Calculate the moments
% Takes 38s for GLAD, 22s for deep GLAD
% Takes 455s for LASER, 88s for deep LASER

tic
for i = 1:length(dist_axis)
    disp(i)
    id = find(dist>= dist_bin(i) & dist<dist_bin(i+1));
    
    pairs_sep(i).dul = dul(id);
    pairs_sep(i).dut = dut(id);
    
    pairs_per_bin(i) = length(id);
    SF1l(i) = nanmean(pairs_sep(i).dul.^1);
    SF1t(i) = nanmean(pairs_sep(i).dut.^1);
        
    SF2ll(i) = nanmean(pairs_sep(i).dul.^2);
    SF2tt(i) = nanmean(pairs_sep(i).dut.^2);
    SF2lt(i) = nanmean(pairs_sep(i).dut.*pairs_sep(i).dul);
    
    SF3lll(i) = nanmean(pairs_sep(i).dul.^3);
    SF3ltt(i) = nanmean(pairs_sep(i).dul.*pairs_sep(i).dut.^2);
end
toc



%% 
clear dist dul dut

save ../data/SF_GLAD_deep_500m.mat SF1l SF1t SF2ll SF2tt SF2lt SF3lll SF3ltt dist_axis pairs_per_bin
%save ../data/SF_LASER_deep_500m.mat SF1l SF1t SF2ll SF2tt SF2lt SF3lll SF3ltt dist_axis pairs_per_bin
%save ../data/SF_LASER_deep_500m_box.mat SF1l SF1t SF2ll SF2tt SF2lt SF3lll SF3ltt dist_axis
%save ../data/SF_GLAD.mat SF2ll SF2tt SF3lll SF3ltt dist_axis
%save ../data/SF_LASER.mat SF2ll SF2tt SF3lll SF3ltt dist_axis

%% Some exploratory plots below
% SF2

figure
loglog(dist_axis, SF2ll, 'linewidth',2)
hold all
loglog(dist_axis, SF2tt, 'linewidth',2)

loglog(dist_axis, SF2ll+SF2tt, 'linewidth',2)


loglog(dist_axis, 1e-4*dist_axis.^(2/3), '--', 'color','k')
loglog(dist_axis, 1e-7*dist_axis.^(2), '--', 'color','k')

%axis([10 1000e3 5e-5 5e-1])

%% SF3 
figure
loglog(dist_axis, abs(SF3lll+SF3ltt), 'linewidth',2)
hold all
loglog(dist_axis, SF3lll+SF3ltt, '+', 'linewidth',2)
loglog(dist_axis, -SF3lll-SF3ltt, 'o', 'linewidth',2)
axis([10 1000e3 1e-8 1])

%% SF3 2 parts
figure
loglog(dist_axis, abs(SF3lll), 'linewidth',2)
hold all
loglog(dist_axis, abs(SF3ltt), 'linewidth',2)
%loglog(dist_axis, SF3lll, '+', 'linewidth',2)
%loglog(dist_axis, -SF3lll, 'o', 'linewidth',2)
axis([10 1000e3 1e-8 1])

%% Compute the histograms with normalized velocities 

norm_vel_bins = linspace(-10,10,51);
dx = norm_vel_bins(2) - norm_vel_bins(1); 
norm_vel_axis = 0.5*( norm_vel_bins(1:end-1) + norm_vel_bins(2:end) );
%
Norm_vel_hist = zeros(length(norm_vel_axis), length(dist_axis)); 
%%
for i=1:length(dist_axis)
   Norm_vel_hist(:, i) = histcounts(pairs_sep(i).dul/ SF2ll(i)^0.5, norm_vel_bins, 'Normalization', 'probability')/dx;
end
   

%% Standard normal 

f_normal = exp(-norm_vel_axis.^2 / 2)/ sqrt(2*pi);
f_normal_rep = repmat(f_normal', 1, 30 );

%%
figure
semilogy(norm_vel_axis, f_normal, '--', 'color','k')
hold all 

semilogy(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=500,1,'last')), 'linewidth',2)
semilogy(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=5000,1,'last')), 'linewidth',2)
semilogy(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=50000,1,'last')), 'linewidth',2)

legend('Gaussian', '500m', '5km', '50km')

xlabel('$\delta u_l/\sigma$', 'interpreter','latex')
ylabel('PDF')
axis([-8 8 5e-4 1])
set(gca,'fontsize',16,'fontname','times')
%print('./figures/GLAD_PDF_zoom.eps','-depsc', '-r400')
print('./figures/LASER_PDF_zoom_box.eps','-depsc', '-r400')
%%
figure
contourf(log10(dist_axis), norm_vel_axis, log10(Norm_vel_hist), 'Linecolor', 'none')
hold on
%grid on
contour(log10(dist_axis), norm_vel_axis, log10(f_normal_rep), [-7, -5, -3, -1], 'linewidth', 3)
contour(log10(dist_axis), norm_vel_axis, log10(f_normal_rep), [-7, -5, -3, -1], '--', 'linecolor', 'k', 'linewidth', 1)

colorbar 
caxis([-7 0])
