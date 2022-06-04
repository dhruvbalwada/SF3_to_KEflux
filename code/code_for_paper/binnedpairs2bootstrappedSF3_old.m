%% Binned pairs to 3rd Order Velocity Structure Functions 
% Dhruv Balwada, 29 January 2021
%
% Modified on 18 May 2021

% This script takes the pairs that are separated in different separation
% bins and estimates the third order structure function 
% SF3 = <du_l (du_l^2 + du_t^2)>
%
% Here do bootstrapping to generate confidence intervals on SF3. The N pairs
% present in each bin are resampled to get N random samples (with
% repetition). These resampled data are then used to create an estimate of
% SF3. This process is repeated many times to generate many SF3s.

%%%%%%%
%%%%%%%
%%%%%%% This does regular bootstrapping
%%%%%%% please use the block bootstrapping code instead.
%%%%%%%
%%%%%%%


clear all
close all
%% Load the file where the binned pairs are
load('../data/structure_pairs_GLAD_deep_500m.mat')
%load('../data/structure_pairs_LASER_deep_500m.mat')

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
%% Calculate the moments
% Takes 38s for GLAD, 22s for deep GLAD
% Takes 455s for LASER, 88s for deep LASER

tic
for i = 1:length(dist_axis)
    disp(i)
    id = find(dist>= dist_bin(i) & dist<dist_bin(i+1));
    
    pairs_sep(i).dul = dul(id);
    pairs_sep(i).dut = dut(id);
    
    %pairs_per_bin(i) = length(id);
    %SF1l(i) = nanmean(pairs_sep(i).dul.^1);
    %SF1t(i) = nanmean(pairs_sep(i).dut.^1);
        
    %SF2ll(i) = nanmean(pairs_sep(i).dul.^2);
    %SF2tt(i) = nanmean(pairs_sep(i).dut.^2);
    %SF2lt(i) = nanmean(pairs_sep(i).dut.*pairs_sep(i).dul);
    
    %SF3lll(i) = nanmean(pairs_sep(i).dul.^3);
    %SF3ltt(i) = nanmean(pairs_sep(i).dul.*pairs_sep(i).dut.^2);
end
toc


%% Generate many estimates of SF3 with random samples

l=length(dist_axis);

bootsamples = 1000; 
s3lll = zeros(l,bootsamples);
s3ltt = zeros(l,bootsamples);

% this loop is extremely extremely slow
% 206 minutes for GLAD
% 
tic
for i =1:l
    
    npairs = length(pairs_sep(i).dul);
    disp(i)
    
    s3lll(i,:) = bootstrp(bootsamples, @mean, pairs_sep(i).dul.^3); 
    s3ltt(i,:) = bootstrp(bootsamples, @mean, pairs_sep(i).dul.*pairs_sep(i).dut.^2 ); 
    %for j =1:bootsamples
    %    y = randsample(npairs, npairs, true);
        
    %    s3lll(i,j) = nanmean(pairs_sep(i).dul(y).^3);
    %    s3ltt(i,j) = nanmean(pairs_sep(i).dul(y).*pairs_sep(i).dut.^2);
    %end
end
toc

%% Save to file

save ../data/GLAD_S3_deep500_boot_strap.mat s3lll s3ltt dist_axis dist_bin
%save GLAD_S3_boot_strap.mat s3lll s3ltt dist_axis dist_bin

%%

load ../data/GLAD_S3_deep500_boot_strap.mat

%%

SF3_mean_bs = mean(s3lll+s3ltt,2); 
SF3_stderr_bs = std(s3lll+s3ltt,0,2); 


