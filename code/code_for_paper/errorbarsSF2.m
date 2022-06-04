%% Estimates of errorbars on SF2
% Dhruv Balwada, 15 December 2021


clear all
%close all
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

%%
figure
loglog(dist_axis, dof)

%%
for i = 1:length(pairs_sep)
    npairs_sep(i) = length(pairs_sep(i).dul);
    n_blocks_sep(i) =  dof(i); % number of blocks at that separation (basically the dof)
    nsamps_per_block_sep(i) = ceil(npairs_sep(i)/ n_blocks_sep(i));
end

%%
if strcmp(experiment, 'GLAD')
    save ../data/GLAD_npairs_deep500_block_boot_strap_Ldof.mat npairs_sep dist_axis dist_bin dof nsamps_per_block_sep
else
    save ../data/LASER_npairs_deep500_box_constrained_block_boot_strap_Ldof.mat npairs_sep dist_axis dist_bin dof nsamps_per_block_sep
end

stop
%%
clear SF3 SF3_mean SF3_stderr

num_boot = 599;
SF2ll = zeros(length(dist_axis), num_boot);
SF2ll_mean = zeros(length(dist_axis),1);
SF2ll_stderr = zeros(length(dist_axis),1);

SF2tt = zeros(length(dist_axis), num_boot);
SF2tt_mean = zeros(length(dist_axis),1);
SF2tt_stderr = zeros(length(dist_axis),1);

%%
tic
for i = 1:length(dist_axis)
    
    
    disp(i)
    blocksize = nsamps_per_block_sep(i);
    %blocksize = 1;
    numblocks = floor(npairs_sep(i)/ blocksize);
    
    if npairs_sep(i)>10
        n = numblocks*blocksize;
        
        blocks_dul = reshape(pairs_sep(i).dul(1:n), [blocksize, numblocks])';
        blocks_dut = reshape(pairs_sep(i).dut(1:n), [blocksize, numblocks])';
        
        SF2ll_samp = blocks_dul.^2; %+ blocks_dul.*blocks_dut.^2;
        SF2tt_samp = blocks_dut.^2;
        
        % create blocks of bootstrap samples
        %SF3_bs = bootstrp(num_boot, @(x)x', SF3_samp');
        % calculate means of each bootstrap sample
        %SF3 = mean(SF3_bs, 2);
        
        [SF2ll(i,:), bootsamp] = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF2ll_samp);
        % the double mean above first takes mean over the blocks, then averages
        % the different blocks.
        % the below loop is needed because we want correspondence between
        % the SF2ll and SF2tt samples
        for j = 1:num_boot
            SF2tt(i,j) = mean(mean(SF2tt_samp(bootsamp(:,j),:),2),1);
        end
        
    else
        SF2ll(i,:) = NaN;
        SF2tt(i,:) = NaN; 
    end
    % Mean and standard error of the estimates
    SF2ll_mean(i)   = mean(SF2ll(i,:));
    SF2ll_stderr(i) = std(SF2ll(i,:)); % boot strap std err is the std of bs estimates
    SF2tt_mean(i)   = mean(SF2tt(i,:));
    SF2tt_stderr(i) = std(SF2tt(i,:));
    
end
toc

%%

if strcmp(experiment, 'GLAD')
    save ../data/GLAD_S2_deep500_block_boot_strap_Ldof.mat SF2ll SF2tt SF2ll_mean SF2ll_stderr SF2tt_mean SF2tt_stderr dof dist_axis dist_bin
else
    save ../data/LASER_S2_deep500_box_constrained_block_boot_strap_Ldof.mat SF2ll SF2tt SF2ll_mean SF2ll_stderr SF2tt_mean SF2tt_stderr dof dist_axis dist_bin
end

%%
ebar(1,:)  =  2*0.434*SF2ll_stderr./SF2ll_mean; 
ebar(2,:)  =  2*0.434*SF2ll_stderr./SF2ll_mean; 

for i = 1:size(SF2ll,1)
    CI_SF2ll(:,i) = prctile(SF2ll(i,:), [95,5]);
    CI_SF2tt(:,i) = prctile(SF2tt(i,:), [95,5]);
end
%%
figure
%shadedErrorBar_log(dist_axis, SF2ll_mean, ebar)
%errorbar(dist_axis, SF2ll_mean./dist_axis', SF2ll_stderr./dist_axis')
shadedErrorBar_log(dist_axis, SF2ll_mean, CI_SF2ll, {'color', 'r'})
hold all
%shadedErrorBar_log(dist_axis, SF2tt_mean, CI_SF2tt, {'color', 'b'})
%shadedErrorBar_log(dist_axis, SF2ll_mean, CI_SF2ll)

%errorbar(dist_axis, SF2tt_mean, SF2tt_stderr/2)
%set(gca, 'XScale','log', 'YScale','linear')

%errorbar(log10(dist_axis), SF3_mean_bs./dist_axis', SF3_stderr_bs./dist_axis'/2)
grid on