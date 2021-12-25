% Dhruv Balwada
% 24 Dec 2021 

% Here we make some plots to show how the 
% modified block bootstrapping works 

clear all
close all 

%%
%GLAD
load('../data/structure_pairs_GLAD_deep_500m.mat');
%LASER = load('../data/structure_pairs_GLAD_deep_500m.mat');

%%
tpts = length(pairs_time);
%
npairs = zeros(tpts,1);
for i = 1:tpts
    npairs(i) = length(find(~isnan(pairs_time(i).dul)));
end
%% % Align pairs in a single vector
% 
dul = zeros(sum(npairs),1);
dut = zeros(sum(npairs),1);
dist = zeros(sum(npairs),1);

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
clear pair_times

%% 
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
% Generate vel axis
vel_bins = linspace(-2, 2, 50);
vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));

%%
edges = linspace(0, 0.1, 101); 
id = 12;
figure
h = histogram(pairs_sep(id).dul.^2, edges, 'Normalization','pdf');
hold all
vlines(mean(pairs_sep(id).dul.^2),'r--')
set(gca, 'Yscale','log')
xlabel('$\delta u_l^2 [m^2/s^2]$', 'interpreter','Latex')
ylabel('PDF')
ylim([.01 1000])
set(gca,'Fontsize', 16)
print(['du2.png'],'-dpng', '-r400')

%%
edges = linspace(-0.01, 0.01, 101); 
id = 12;
figure
h = histogram(pairs_sep(id).dul.^3, edges, 'Normalization','pdf');
hold all
vlines(mean(pairs_sep(id).dul.^3),'r--')
set(gca, 'Yscale','log')
xlabel('$\delta u_l^3 [m^3/s^3]$', 'interpreter','Latex')
ylabel('PDF')
ylim([.01 4000])
set(gca,'Fontsize', 16)
print(['du3.png'],'-dpng', '-r400')

%%
id = 12;
figure
plot(pairs_sep(id).dul(1:500))
xlabel('Sample number')
set(gca,'Fontsize', 16,'Fontname','Times')

ylabel('$\delta u (r)$','interpreter','Latex')
title(['Bin:', num2str(ceil(dist_bin(id))), '-', num2str(ceil(dist_bin(id+1))), 'm'])
print(['samples1.png'],'-dpng', '-r400')

%%
id = 21;
figure
plot(pairs_sep(id).dul(1:500))
xlabel('Sample number')
set(gca,'Fontsize', 16,'Fontname','Times')

ylabel('$\delta u (r)$','interpreter','Latex')
title(['Bin:', num2str(ceil(dist_bin(id))), '-', num2str(ceil(dist_bin(id+1))), 'm'])
print(['samples2.png'],'-dpng', '-r400')

%%

GLAD_dof = load('../data/GLAD_npairs_deep500_block_boot_strap_Ldof.mat');
LASER_dof = load('../data/LASER_npairs_deep500_box_constrained_block_boot_strap_Ldof.mat');


%%
cols=colororder;
%%

figure
loglog(GLAD_dof.dist_axis/1e3, GLAD_dof.dof,...
    'linewidth',2, 'color',cols(7,:))
hold all
loglog(LASER_dof.dist_axis/1e3, LASER_dof.dof, ...
    'linewidth',2, 'color',cols(1,:))

axis([0.01 1000 1 1e4])
set(gca,'fontsize',16, 'fontname','Times')
legend( 'GLAD', 'LASER', 'location','northeast')

ylabel('$N^{DOF}$', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')

print(['ndof.png'],'-dpng', '-r400')

%%
figure
loglog(GLAD_dof.dist_axis/1e3, GLAD_dof.nsamps_per_block_sep,...
    'linewidth',2, 'color',cols(7,:))
hold all
loglog(LASER_dof.dist_axis/1e3, LASER_dof.nsamps_per_block_sep, ...
    'linewidth',2, 'color',cols(1,:))

axis([0.01 1000 1 1e7])
set(gca,'fontsize',16, 'fontname','Times')
legend( 'GLAD', 'LASER', 'location','northwest')

ylabel('Number of samples per block', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')

print(['nsamps_per_block.png'],'-dpng', '-r400')