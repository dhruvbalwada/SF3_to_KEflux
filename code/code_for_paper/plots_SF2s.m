% Plots of SF2 with errorbars
% Dhruv Balwada
% 19 December 2021

clear all
close all 

%% Load data

GLAD = load('../data/GLAD_S2_deep500_block_boot_strap_Ldof.mat');
LASER = load('../data/LASER_S2_deep500_box_constrained_block_boot_strap_Ldof.mat');

%%
GLAD.SF2tot  = GLAD.SF2ll  + GLAD.SF2tt; 
LASER.SF2tot = LASER.SF2ll + LASER.SF2tt; 

for i = 1:size(GLAD.SF2tot,1)
    GLAD.CI_SF2tot(:,i) = prctile(GLAD.SF2tot(i,:), [95,5]);
    LASER.CI_SF2tot(:,i) = prctile(LASER.SF2tot(i,:), [95,5]);
    
    GLAD.CI_SF2tt(:,i) = prctile(GLAD.SF2tt(i,:), [95,5]);
    LASER.CI_SF2tt(:,i) = prctile(LASER.SF2tt(i,:), [95,5]);
    
    GLAD.CI_SF2ll(:,i) = prctile(GLAD.SF2ll(i,:), [95,5]);
    LASER.CI_SF2ll(:,i) = prctile(LASER.SF2ll(i,:), [95,5]);
    
end

GLAD.SF2tot_mean = mean(GLAD.SF2tot, 2);
LASER.SF2tot_mean = mean(LASER.SF2tot, 2);

GLAD.SF2tt_mean = mean(GLAD.SF2tt, 2);
LASER.SF2tt_mean = mean(LASER.SF2tt, 2);

GLAD.SF2ll_mean = mean(GLAD.SF2ll, 2);
LASER.SF2ll_mean = mean(LASER.SF2ll, 2);

%%
cols = colororder;
r = GLAD.dist_axis/1e3;
%% Total SF2 plot
figure
G = shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2tot_mean, ...
    GLAD.CI_SF2tot, {'color', cols(7,:), 'linewidth',2}, 1);
hold all
L = shadedErrorBar_log(LASER.dist_axis/1e3, LASER.SF2tot_mean, ...
    LASER.CI_SF2tot, {'color', cols(1,:), 'linewidth',2}, 1);

r2 = loglog(r, 5e-3*r.^2, 'k--');
r23 = loglog(r, 5e-3*r.^(2/3), 'k-');

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',18, 'fontname','Times')
%legend([G.mainLine, L.mainLine, r2, r23], '{\color[rgb]{0.635,0.078,0.184}} GLAD/ Summer', ...
%    '\color[rgb]{0,0.447,0.741} LASER/ Winter', '$r^{2}$', '$r^{2/3}$', 'location','northwest', 'interpreter','latex')
legend([G.mainLine, L.mainLine, r2, r23], 'GLAD/ Summer', 'LASER/ Winter', '$r^{2}$', '$r^{2/3}$', 'location','northwest', 'interpreter','latex')

ylabel('$<\delta u_{ll}^2 + \delta u_{tt}^2> [m^2/s^2]$', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')

print(['SF2_tot_both.png'],'-dpng', '-r400')




%% Ro plot
f = 2 * (2*pi/24/3600) * sind(27); 

GLAD.Ro = sqrt(GLAD.SF2tot_mean)/f./GLAD.dist_axis'; 
LASER.Ro = sqrt(LASER.SF2tot_mean)/f./LASER.dist_axis'; 

GLAD.CI_Ro = sqrt(GLAD.CI_SF2tot)/f./GLAD.dist_axis; 
LASER.CI_Ro = sqrt(LASER.CI_SF2tot)/f./LASER.dist_axis; 

figure
G = shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.Ro, ...
    GLAD.CI_Ro, {'color', cols(7,:), 'linewidth',2}, 1);
hold all
L = shadedErrorBar_log(LASER.dist_axis/1e3, LASER.Ro, ...
    LASER.CI_Ro, {'color', cols(1,:), 'linewidth',2}, 1);
hlines(1, '1--')
axis([0.01 1000 0.01 20])

set(gca,'fontsize',20, 'fontname','Times')
%legend('GLAD', 'LASER','Ro=1', 'location', 'best')
ylabel('$Ro $', 'interpreter', 'latex')
xlabel('$r [km]$', 'interpreter', 'latex')

print(['rossby_num.png'],'-dpng', '-r400')

%% SF2  GLAD components

figure
R = shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2ll_mean, GLAD.CI_SF2ll,...
    {'color',cols(1,:), 'linewidth',2},1);
hold all


D = shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2tt_mean, GLAD.CI_SF2tt,...
    {'color',cols(2,:), 'linewidth',2}, 1);


E = shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2tot_mean, GLAD.CI_SF2tot,...
    {'color',cols(3,:), 'linewidth',2}, 1);
%shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2tot_mean, GLAD.CI_SF2tot,...
%    {'color',cols(1,:), 'linewidth',2}, 2)

r2 = loglog(r, 5e-3*r.^2, 'k--');
r23 = loglog(r, 5e-3*r.^(2/3), 'k-');

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',18, 'fontname','Times')
legend([R.mainLine, D.mainLine, E.mainLine, r2, r23], 'Longitudinal', 'Transverse', 'Total',  '$r^{2}$', '$r^{2/3}$', ...
    'location','northwest', 'interpreter', 'latex')
%legend('GLAD rot', 'GLAD div', 'GLAD total', 'LASER rot', 'LASER div', 'LASER total' , 'location', 'best')
ylabel('$<\delta u_i^2> [m^2/s^2]$', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')
title('GLAD/ Summer', 'color', cols(7,:))
print(['GLAD_SF2_components.png'],'-dpng', '-r400')

%% SF2 LASER components

figure
R = shadedErrorBar_log(LASER.dist_axis/1e3, LASER.SF2ll_mean, LASER.CI_SF2ll,...
    {'color',cols(1,:), 'linewidth',2},1);
hold all


D = shadedErrorBar_log(LASER.dist_axis/1e3, LASER.SF2tt_mean, LASER.CI_SF2tt,...
    {'color',cols(2,:), 'linewidth',2}, 1);


E = shadedErrorBar_log(LASER.dist_axis/1e3, LASER.SF2tot_mean, LASER.CI_SF2tot,...
    {'color',cols(3,:), 'linewidth',2}, 1);
%shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2tot_mean, GLAD.CI_SF2tot,...
%    {'color',cols(1,:), 'linewidth',2}, 2)

r2 = loglog(r, 5e-3*r.^2, 'k--');
r23 = loglog(r, 5e-3*r.^(2/3), 'k-');

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',18, 'fontname','Times')
legend([R.mainLine, D.mainLine, E.mainLine, r2, r23], 'Longitudinal', 'Transverse', 'Total',  '$r^{2}$', '$r^{2/3}$', ...
    'location','northwest', 'interpreter', 'latex')
%legend('GLAD rot', 'GLAD div', 'GLAD total', 'LASER rot', 'LASER div', 'LASER total' , 'location', 'best')
ylabel('$<\delta u_i^2> [m^2/s^2]$', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')
title('LASER/ Winter', 'color', cols(1,:))
print(['LASER_SF2_components.png'],'-dpng', '-r400')

%% Helmholtz decomposition 

GLAD.SF2rr = 0*GLAD.SF2ll;
GLAD.SF2dd = 0*GLAD.SF2ll;

LASER.SF2rr = 0*GLAD.SF2ll;
LASER.SF2dd = 0*GLAD.SF2ll;

for i=1:size(GLAD.SF2ll,2) 
    [GLAD.SF2rr(:,i), GLAD.SF2dd(:,i)] = helmholtz_decompose(GLAD.dist_axis', GLAD.SF2ll(:,i), GLAD.SF2tt(:,i)); 
    [LASER.SF2rr(:,i), LASER.SF2dd(:,i)] = helmholtz_decompose(LASER.dist_axis', LASER.SF2ll(:,i), LASER.SF2tt(:,i));
end


for i = 1:size(GLAD.SF2rr,1)
    GLAD.CI_SF2rr(:,i)  = prctile(GLAD.SF2rr(i,:), [95,5]);
    LASER.CI_SF2rr(:,i) = prctile(LASER.SF2rr(i,:), [95,5]);
    GLAD.CI_SF2dd(:,i)  = prctile(GLAD.SF2dd(i,:), [95,5]);
    LASER.CI_SF2dd(:,i) = prctile(LASER.SF2dd(i,:), [95,5]);
end

GLAD.SF2rr_mean = mean(GLAD.SF2rr, 2);
LASER.SF2rr_mean = mean(LASER.SF2rr, 2);
GLAD.SF2dd_mean = mean(GLAD.SF2dd, 2);
LASER.SF2dd_mean = mean(LASER.SF2dd, 2);
disp('done')


%% SF2 Helmholtz GLAD

figure
R = shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2rr_mean, GLAD.CI_SF2rr,...
    {'color',cols(4,:), 'linewidth',2},1);
hold all

GLAD.SF2dd_mean(GLAD.SF2dd_mean<=0) = 1e-7;
GLAD.CI_SF2dd(GLAD.CI_SF2dd<=0) = 1e-7;

D = shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2dd_mean, GLAD.CI_SF2dd,...
    {'color',cols(3,:), 'linewidth',2}, 1);
%shadedErrorBar_log(GLAD.dist_axis/1e3, GLAD.SF2tot_mean, GLAD.CI_SF2tot,...
%    {'color',cols(1,:), 'linewidth',2}, 2)
r2 = loglog(r, 5e-3*r.^2, 'k--');
r23 = loglog(r, 5e-3*r.^(2/3), 'k-');

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',18, 'fontname','Times')
legend([R.mainLine, D.mainLine], 'Rotational', 'Divergent', ...
    'location','northwest')
%legend('GLAD rot', 'GLAD div', 'GLAD total', 'LASER rot', 'LASER div', 'LASER total' , 'location', 'best')
ylabel('$<\delta u_i^2> [m^2/s^2]$', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')

print(['GLAD_SF2_decompose.png'],'-dpng', '-r400')

%% LASER
figure
shadedErrorBar_log(LASER.dist_axis/1e3, LASER.SF2rr_mean, LASER.CI_SF2rr,...
    {'color',cols(4,:), 'linewidth',2},1);
hold all

LASER.SF2dd_mean(LASER.SF2dd_mean<=0) = 1e-7;
LASER.CI_SF2dd(LASER.CI_SF2dd<=0) = 1e-7;

shadedErrorBar_log(LASER.dist_axis/1e3, LASER.SF2dd_mean, LASER.CI_SF2dd,...
    {'color',cols(3,:), 'linewidth',2}, 1);
r2 = loglog(r, 5e-3*r.^2, 'k--');
r23 = loglog(r, 5e-3*r.^(2/3), 'k-');
%shadedErrorBar_log(LASER.dist_axis/1e3, LASER.SF2tot_mean, LASER.CI_SF2tot,...
%    {'color',cols(1,:), 'linewidth',2}, 2)

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',18, 'fontname','Times')
%legend('GLAD rot', 'GLAD div', 'GLAD total', 'LASER rot', 'LASER div', 'LASER total' , 'location', 'best')
ylabel('$<\delta u_i^2> [m^2/s^2]$', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')

print(['LASER_SF2_decompose.png'],'-dpng', '-r400')

%% Ratios 

Ratio_tot = LASER.SF2tot./GLAD.SF2tot; 
Ratio_rr  = LASER.SF2rr./GLAD.SF2rr; 

Ratio_tot_mean = mean(Ratio_tot, 2); 
Ratio_rr_mean = mean(Ratio_rr, 2); 

for i = 1:size(Ratio_tot,1)
    CI_Ratio_tot(:,i) = prctile(Ratio_tot(i,:), [95,5]);
    CI_Ratio_rr(:,i)  = prctile(Ratio_rr(i,:), [95,5]);
end

%% Plot of ratios 

figure
T = shadedErrorBar_semilogx(LASER.dist_axis/1e3, Ratio_tot_mean, CI_Ratio_tot,...
    {'color',cols(5,:), 'linewidth',2},1);
hold all
R = shadedErrorBar_semilogx(LASER.dist_axis/1e3, Ratio_rr_mean, CI_Ratio_rr,...
    {'color',cols(4,:), 'linewidth',2},1);

hlines(1, '1--')
axis([0.01 1000 0.01 18])
yticks([1,6,11,17])
legend([T.mainLine, R.mainLine], 'Total', 'Rotational')
set(gca,'fontsize',18, 'fontname','Times')
xlabel('$r [km]$', 'interpreter','latex')
ylabel('LASER/GLAD')

print(['ratio.png'],'-dpng', '-r400')

 %% GLAD componensated decompose
figure
R = shadedErrorBar_semilogx(GLAD.dist_axis/1e3, GLAD.SF2rr_mean./GLAD.dist_axis', GLAD.CI_SF2rr./GLAD.dist_axis,...
    {'color',cols(4,:), 'linewidth',2},1);
hold all

D = shadedErrorBar_semilogx(GLAD.dist_axis/1e3, GLAD.SF2dd_mean./GLAD.dist_axis', GLAD.CI_SF2dd./GLAD.dist_axis,...
    {'color',cols(3,:), 'linewidth',2},1);

T = shadedErrorBar_semilogx(GLAD.dist_axis/1e3, GLAD.SF2tot_mean./GLAD.dist_axis', GLAD.CI_SF2tot./GLAD.dist_axis,...
    {'color',cols(5,:), 'linewidth',2},1);
hlines(0, '--')
grid on
axis([0.01 1000 -3e-6 14e-6])
set(gca,'fontsize',16, 'fontname','Times')
legend([R.mainLine, D.mainLine, T.mainLine], 'Rotational', 'Divergent', 'Total', 'location', 'best')
ylabel('$<\delta u_i^2>/r [m/s^2]$', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')
title('GLAD/ Summer', 'color', cols(7,:))

print(['GLAD_SF2_decompose_compensated.png'],'-dpng', '-r400')

 %% LASER componensated decompose
figure
R = shadedErrorBar_semilogx(LASER.dist_axis/1e3, LASER.SF2rr_mean./LASER.dist_axis', LASER.CI_SF2rr./LASER.dist_axis,...
    {'color',cols(4,:), 'linewidth',2},1);
hold all

D = shadedErrorBar_semilogx(LASER.dist_axis/1e3, LASER.SF2dd_mean./LASER.dist_axis', LASER.CI_SF2dd./LASER.dist_axis,...
    {'color',cols(3,:), 'linewidth',2},1);

T = shadedErrorBar_semilogx(LASER.dist_axis/1e3, LASER.SF2tot_mean./LASER.dist_axis', LASER.CI_SF2tot./LASER.dist_axis,...
    {'color',cols(5,:), 'linewidth',2},1);

%shadedErrorBar_log(LASER.dist_axis/1e3, LASER.SF2tot_mean, LASER.CI_SF2tot,...
%    {'color',cols(1,:), 'linewidth',2}, 2)
hlines(0, '--')
grid on
axis([0.01 1000 -3e-6 14e-6])
set(gca,'fontsize',16, 'fontname','Times')
%legend([R.mainLine, D.mainLine, T.mainLine], 'Rotational', 'Divergent', 'Total', 'location', 'best')
ylabel('$<\delta u_i^2>/r [m/s^2]$', 'interpreter','latex')
xlabel('$r [km]$', 'interpreter','latex')
title('LASER/ Winter', 'color', cols(1,:))

print(['LASER_SF2_decompose_compensated.png'],'-dpng', '-r400')

