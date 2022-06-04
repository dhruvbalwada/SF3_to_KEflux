%% Estimates of errorbars on SF3
% Dhruv Balwada, 15 December 2021

% The primary purpose of this script is to try and figure
% out how to estimate reasonable errorbars on SF3.

% The problem with standard bootstrapping is that the
% samples are not independent - the number of DOF
% are much smaller than the number of pairs, because
% samples are correlated in time.

clear all
close all
%% Load the file where the binned pairs are

experiment = 'LASER';

if strcmp(experiment, 'GLAD')
    load('../data/structure_pairs_GLAD_deep_500m.mat')
    Ttot = 90*24*3600;
else
    load('../data/structure_pairs_LASER_deep_500m_box_constrained.mat')
    Ttot = 60*24*3600; % shorter exp duration
end

%% Get time axis and figure out how many pairs
tpts = length(pairs_time);
%
npairs = zeros(tpts,1);
for i = 1:tpts
    npairs(i) = length(find(~isnan(pairs_time(i).dul)));
end

%% Make into a single vector
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
%%
tic
for i = 1:length(dist_axis)
    disp(i)
    id = find(dist>= dist_bin(i) & dist<dist_bin(i+1));
    
    pairs_sep(i).dul = dul(id);
    pairs_sep(i).dut = dut(id);
end
toc
%%
% Compute mean SF2 to use for estimating DOF

for i = 1:length(dist_axis)
    %pairs_per_bin(i) = length(id);
    %SF1l(i) = nanmean(pairs_sep(i).dul.^1);
    %SF1t(i) = nanmean(pairs_sep(i).dut.^1);
    
    SF2ll(i) = nanmean(pairs_sep(i).dul.^2);
    SF2tt(i) = nanmean(pairs_sep(i).dut.^2);
    %SF2lt(i) = nanmean(pairs_sep(i).dut.*pairs_sep(i).dul);
    
    %SF3lll(i) = nanmean(pairs_sep(i).dul.^3);
    %SF3ltt(i) = nanmean(pairs_sep(i).dul.*pairs_sep(i).dut.^2);
end


%% test for setting up block bootstrap

test_flag =0 ;
if test_flag == 1
    ts = 1:12;
    blockSize = 2;
    numBlocks = length(ts) / blockSize;           % must be integer
    %blocks = reshape(ts, [numBlocks,blockSize])  % reshape into non-overlapping blocks
    blocks = reshape(ts, [blockSize, numBlocks])';
    nSamples = 10;
    samples = bootstrp(1, @(x)x', blocks);
    % the funny x' thing happens because the data is being converted to a row vector
    
end
%% Degree of freedom using time of process and total length of experiment
%
%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
Tscale_tot = 1./(((SF2ll +SF2tt).^0.5)./dist_axis);
Tscale_ll = 1./(((SF2ll).^0.5)./dist_axis);
Tscale_tt = 1./(((SF2tt).^0.5)./dist_axis);

dof = ceil(Ttot./Tscale_tot); % this is essentially T_tot/T_scale(r)

%% Plot of time scales
figure
loglog(dist_axis/1e3, Tscale_ll/24/3600)
hold all
loglog(dist_axis/1e3, Tscale_tt/24/3600)
loglog(dist_axis/1e3, Tscale_tot/24/3600)

%% Plot of degrees of freedom (contained in the total time)
figure
loglog(dist_axis/1e3, dof)

%%
for i = 1:length(pairs_sep)
    npairs_sep(i) = length(pairs_sep(i).dul);
    n_blocks_sep(i) =  dof(i); % number of blocks at that separation (basically the dof)
    nsamps_per_block_sep(i) = ceil(npairs_sep(i)/ n_blocks_sep(i));
end

%%
figure,
loglog(dist_axis, nsamps_per_block_sep)
hold all
loglog(dist_axis, npairs_sep)
legend('Pairs per block', 'Pairs')

%%
% ~450 seconds for GLAD

clear SF3 SF3_mean SF3_stderr

num_boot = 599;
SF3 = zeros(length(dist_axis), num_boot);
SF3_mean = zeros(length(dist_axis),1);
SF3_stderr = zeros(length(dist_axis),1);

%%
tic
for i = 1:length(dist_axis)
    
    
    disp(i)
    blocksize = nsamps_per_block_sep(i);
    %blocksize = 1;
    numblocks = floor(npairs_sep(i)/ blocksize);
    
    if npairs_sep>10
        n = numblocks*blocksize;
        
        blocks_dul = reshape(pairs_sep(i).dul(1:n), [blocksize, numblocks])';
        blocks_dut = reshape(pairs_sep(i).dut(1:n), [blocksize, numblocks])';
        
        SF3_samp = blocks_dul.^3 + blocks_dul.*blocks_dut.^2;
        
        % create blocks of bootstrap samples
        %SF3_bs = bootstrp(num_boot, @(x)x', SF3_samp');
        % calculate means of each bootstrap sample
        %SF3 = mean(SF3_bs, 2);
        
        SF3(i,:) = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF3_samp);
        % the double mean above first takes mean over the blocks, then averages
        % the different blocks.
    else
        SF3(i,:) = NaN;
    end
    % Mean and standard error of the estimates
    SF3_mean(i) = mean(SF3(i,:));
    SF3_stderr(i) = std(SF3(i,:)); % boot strap std err is the std of bs estimates
    
end
toc

%%

if strcmp(experiment, 'GLAD')
    save ../data/GLAD_S3_deep500_block_boot_strap_Ldof.mat SF3 SF3_mean SF3_stderr dof dist_axis dist_bin
else
    save ../data/LASER_S3_deep500_box_constrained_block_boot_strap_Ldof.mat SF3 SF3_mean SF3_stderr dof dist_axis dist_bin
end


%%
load ../data/LASER_S3_deep500_box_constrained_block_boot_strap_Ldof.mat
figure
errorbar(log10(dist_axis), SF3_mean./dist_axis', SF3_stderr./dist_axis'/2)
hold all
%errorbar(log10(dist_axis), SF3_mean_bs./dist_axis', SF3_stderr_bs./dist_axis'/2)
grid on

%%
stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Old code below %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Stay out %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Degree of freedom using time of process and total length of time series
% This method does not work well, as it does not account for the fact that
% a lot of timeseries are spatially correlated (particularly at larger
% scales).
%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


%% Time scale (hours)
% Calc a rought time scale using the SF2 (T = r/sqrt(SF2))
Tscalell = ((dist_axis.^2./(SF2ll)).^0.5)/3600;
Tscalett = ((dist_axis.^2./(SF2tt)).^0.5)/3600;
Tscaletot = ((dist_axis.^2./((SF2ll + SF2tt)/2)).^0.5)/3600;

%%
figure
loglog(dist_axis/1e3, Tscalell/24)
hold all
loglog(dist_axis/1e3, Tscalett/24)
loglog(dist_axis/1e3, Tscaletot/24)
grid on
legend('lon', 'tran', 'Tot')

%% Block bootstrap

block_size = ceil(Tscaletot* 4);
% factor of 4 corresponds to the time scale in the data of 15mins
% divide by 2 because SF2 = 2*u^2

%%
for i = 1:length(pairs_sep)
    npairs_sep(i) = length(pairs_sep(i).dul);
    nsamps_sep(i) = npairs_sep(i) / block_size(i);
end

%%
figure,
loglog(dist_axis, npairs_sep)
hold all
loglog(dist_axis, nsamps_sep)
legend('Pairs', 'Samples')

%%
% 568.364847 seconds for GLAD
clear SF3_mean SF3_stderr
SF3_mean = zeros(length(dist_axis),1);
SF3_stderr = zeros(length(dist_axis),1);
num_boot = 599;
%%
tic
for i = 1:length(dist_axis)
    
    disp(i)
    blocksize = block_size(i);
    %blocksize = 1;
    numblocks = floor(npairs_sep(i)/ blocksize);
    
    n = numblocks*blocksize;
    
    blocks_dul = reshape(pairs_sep(i).dul(1:n), [blocksize, numblocks])';
    blocks_dut = reshape(pairs_sep(i).dut(1:n), [blocksize, numblocks])';
    
    SF3_samp = blocks_dul.^3 + blocks_dul.*blocks_dut.^2;
    
    % create blocks of bootstrap samples
    %SF3_bs = bootstrp(num_boot, @(x)x', SF3_samp');
    % calculate means of each bootstrap sample
    %SF3 = mean(SF3_bs, 2);
    
    SF3 = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF3_samp);
    % the double mean above first takes mean over the blocks, then averages
    % the different blocks.
    
    % Mean and standard error of the estimates
    SF3_mean(i) = mean(SF3);
    SF3_stderr(i) = std(SF3); % boot strap std err is the std of bs estimates
    
end
toc

save ../data/GLAD_S3_deep500_block_boot_strap.mat SF3_mean SF3_stderr dist_axis dist_bin

%%
% Even the errorbars made above are not quite right, and likely
% underestimates of true error. This is because there is spatial
% correlation between pairs, which has not been considered. Different parts
% of the timeseries are correlated because they are spatially together,
% while we have only tried to remove the temporal correlations.

%%
% The errorbars have become larger than if we take every sample as
% independent
figure
errorbar(log10(dist_axis), SF3_mean./dist_axis', SF3_stderr./dist_axis'/2)
grid on


%%


%% Older Scraps
%% Compute the histograms with normalized velocities
% Here we tried to calculate what the dist of du^3 look like, hoping to get
% some clues on how to calculate standard errors.
% It turned out
% that they were very long tailed, and not like any standard distributions
% that I could find.

norm_vel_bins = linspace(-0.001,0.001,51);
dx = norm_vel_bins(2) - norm_vel_bins(1);
norm_vel_axis = 0.5*( norm_vel_bins(1:end-1) + norm_vel_bins(2:end) );
%
Norm_vel_hist = zeros(length(norm_vel_axis), length(dist_axis));
%%
for i=1:length(dist_axis)
    Norm_vel_hist(:, i) = histcounts(pairs_sep(i).dul.^3, norm_vel_bins, 'Normalization', 'probability')/dx;
end


%%
figure
semilogy(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=500,1,'last')), 'linewidth',2)
%semilogy(norm_vel_axis, f_normal, '--', 'color','k')
hold all

semilogy(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=5000,1,'last')), 'linewidth',2)
semilogy(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=50000,1,'last')), 'linewidth',2)

legend( '500m', '5km', '50km')

xlabel('$\delta u_l^3$', 'interpreter','latex')
ylabel('PDF')
%axis([-8 8 5e-4 1])
set(gca,'fontsize',16,'fontname','times')

%%
figure
plot(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=500,1,'last')), 'linewidth',2)
%semilogy(norm_vel_axis, f_normal, '--', 'color','k')
hold all

plot(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=5000,1,'last')), 'linewidth',2)
plot(norm_vel_axis, Norm_vel_hist(:,find(dist_axis<=50000,1,'last')), 'linewidth',2)

legend( '500m', '5km', '50km')

xlabel('$\delta u_l^3$', 'interpreter','latex')
ylabel('PDF')
%axis([-8 8 5e-4 1])
set(gca,'fontsize',16,'fontname','times')