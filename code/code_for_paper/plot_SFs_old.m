% Dhruv Balwada
% 16 May 2021

% Script make the plots from both the

clear all
close all

%% Load data

GLAD = load('../data/SF_GLAD_deep_500m.mat');
LASER = load('../data/SF_LASER_deep_500m.mat');
%LASER = load('../data/SF_LASER_deep_500m_box.mat');	

%%

cols = colororder;

%% Pairs per bin

figure 
loglog(GLAD.dist_axis/1e3, GLAD.pairs_per_bin, 'color',cols(1,:), 'linewidth',2)
hold all 
loglog(LASER.dist_axis/1e3, LASER.pairs_per_bin, 'color',cols(2,:), 'linewidth',2)

xlabel('r [km]')
ylabel('Num of pairs')

legend('GLAD', 'LASER', 'location', 'best')

axis([0.01 1000 5e3  1e8])

set(gca,'fontsize',16, 'fontname','Times')
print('./figures/bins_pairs.eps','-depsc', '-r400')
%% SF1

figure
semilogx(GLAD.dist_axis/1e3, GLAD.SF1l, 'color',cols(1,:), 'linewidth',2)
hold all
semilogx(GLAD.dist_axis/1e3, GLAD.SF1t, 'color',cols(2,:), 'linewidth',2)
semilogx(GLAD.dist_axis/1e3, LASER.SF1l, '--', 'color',cols(1,:), 'linewidth',2)
semilogx(GLAD.dist_axis/1e3, LASER.SF1t, '--', 'color',cols(2,:), 'linewidth',2)
xlabel('r [km]')
ylabel('<\delta u_i>')
axis([0.01 1000 -0.01 0.2])
legend('GLAD <\delta u_l>', 'GLAD <\delta u_t>', 'LASER <\delta u_l>', 'LASER <\delta u_t>','location', 'best')
set(gca,'fontsize',16, 'fontname','Times')
%print('./figures/SF1.eps','-depsc', '-r400')


%% SF2 GLAD
figure
loglog(GLAD.dist_axis/1e3, GLAD.SF2ll, 'color',cols(1,:), 'linewidth',2)
hold all
loglog(GLAD.dist_axis/1e3, GLAD.SF2tt, 'color',cols(2,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, GLAD.SF2ll + GLAD.SF2tt, 'color',cols(3,:), 'linewidth',2)


ao = 10^-2;
a1 = 10^-2;
xax = [0.001:10:10^3];
yax2 = ao*(xax.^2);
yax1 = 1.3*a1*(xax.^1);
yax23 = 1.4*a1*(xax.^(2/3));
loglog(xax, yax2,'--','color',[0.5 0.5 0.5],'linewidth',0.7)
%loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',0.7)
loglog(xax, yax23,'-','color',[0.5 0.5 0.5],'linewidth',0.7)
%axis([10^-2 10^3 10^-4 6*10^-2])

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',16, 'fontname','Times')
legend('Longitudinal', 'Transverse','Total', 'r^2', 'r^{2/3}', 'location', 'best')
ylabel('<\delta u_i^2> [m^2/s^2]')
xlabel('r [km]')
print('./figures/SF2_GLAD.eps','-depsc', '-r400')

%% SF2 LASER

figure
loglog(GLAD.dist_axis/1e3, LASER.SF2ll, 'color',cols(1,:), 'linewidth',2)
hold all
loglog(GLAD.dist_axis/1e3, LASER.SF2tt, 'color',cols(2,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, LASER.SF2ll + LASER.SF2tt, 'color',cols(3,:), 'linewidth',2)

loglog(GLAD.dist_axis/1e3, GLAD.SF2ll + GLAD.SF2tt, '--', 'color',cols(3,:), 'linewidth',0.8)

ao = 10^-2;
a1 = 10^-2;
xax = [0.001:10:10^3];
yax2 = ao*(xax.^2);
yax1 = 1.3*a1*(xax.^1);
yax23 = 1.4*a1*(xax.^(2/3));
loglog(xax, yax2,'--','color',[0.5 0.5 0.5],'linewidth',0.7)
%loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',0.7)
loglog(xax, yax23,'-','color',[0.5 0.5 0.5],'linewidth',0.7)
%axis([10^-2 10^3 10^-4 6*10^-2])

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',16, 'fontname','Times')
legend('Longitudinal', 'Transverse','Total', 'GLAD Total', 'r^2', 'r^{2/3}', 'location', 'best')
ylabel('<\delta u_i^2> [m^2/s^2]')
xlabel('r [km]')
print('./figures/SF2_LASER_box.eps','-depsc', '-r400')

%% Ro GLAD and LASER 
f = 2 * (2*pi/24/3600) * sind(27); 
Ro_GLAD = sqrt(GLAD.SF2ll + GLAD.SF2tt)/f./GLAD.dist_axis;
Ro_LASER = sqrt(LASER.SF2ll + LASER.SF2tt)/f./LASER.dist_axis;

figure
loglog(GLAD.dist_axis/1e3, Ro_GLAD, 'color',cols(4,:), 'linewidth',2)
hold all
loglog(LASER.dist_axis/1e3, Ro_LASER, 'color',cols(5,:), 'linewidth',2)
hlines(1, '1b--')
axis([0.01 1000 0.01 20])


set(gca,'fontsize',16, 'fontname','Times')
legend('GLAD', 'LASER','Ro=1', 'location', 'best')
ylabel('$Ro = \sqrt{D2_{ll} + D2_{tt}}/fr$', 'interpreter', 'latex')
xlabel('r [km]')

print('./figures/Ro_box.eps','-depsc', '-r400')

%% Ratios D2 
figure
semilogx(GLAD.dist_axis/1e3, GLAD.SF2tt./ GLAD.SF2ll, 'color',cols(4,:), 'linewidth',2)
hold all
semilogx(LASER.dist_axis/1e3, LASER.SF2tt./ LASER.SF2ll, 'color',cols(5,:), 'linewidth',2)
hlines([3/5, 3], '1b--')
axis([0.01 1000  0, 3.5])


set(gca,'fontsize',16, 'fontname','Times')
legend('GLAD', 'LASER','3/5 - "Divergent+Wave"', '3 - "Rotational+Geostrophic"', 'location', 'best')
ylabel('$Ro = \sqrt{D2_{ll} + D2_{tt}}/fr$', 'interpreter', 'latex')
xlabel('r [km]')

print('./figures/Ratios_box.eps','-depsc', '-r400')

%
%% SF2 Corrected
GLAD.SF2llcor = GLAD.SF2ll - GLAD.SF1l.^2;
GLAD.SF2ttcor = GLAD.SF2tt - GLAD.SF1t.^2;

LASER.SF2llcor = LASER.SF2ll - LASER.SF1l.^2;
LASER.SF2ttcor = LASER.SF2tt - LASER.SF1t.^2;

%% Helmholtz Decomposition

[GLAD.SF2rr, GLAD.SF2dd] = helmholtz_decompose(GLAD.dist_axis, GLAD.SF2ll, GLAD.SF2tt); 
[GLAD.SF2rrcor, GLAD.SF2ddcor] = helmholtz_decompose(GLAD.dist_axis, GLAD.SF2llcor, GLAD.SF2ttcor); 
[LASER.SF2rr, LASER.SF2dd] = helmholtz_decompose(LASER.dist_axis, LASER.SF2ll, LASER.SF2tt); 
[LASER.SF2rrcor, LASER.SF2ddcor] = helmholtz_decompose(LASER.dist_axis, LASER.SF2llcor, LASER.SF2ttcor); 


%% SF2 Helm GLAD
figure
loglog(GLAD.dist_axis/1e3, GLAD.SF2rr, 'color',cols(4,:), 'linewidth',2)
hold all
loglog(GLAD.dist_axis/1e3, GLAD.SF2dd, 'color',cols(5,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, -GLAD.SF2dd, '--','color',cols(5,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, GLAD.SF2rr + GLAD.SF2dd, 'color',cols(3,:), 'linewidth',2)


ao = 10^-2;
a1 = 10^-2;
xax = [0.001:10:10^3];
yax2 = ao*(xax.^2);
yax1 = 1.3*a1*(xax.^1);
yax23 = 1.4*a1*(xax.^(2/3));
loglog(xax, yax2,'--','color',[0.5 0.5 0.5],'linewidth',0.7)
%loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',0.7)
loglog(xax, yax23,'-','color',[0.5 0.5 0.5],'linewidth',0.7)
%axis([10^-2 10^3 10^-4 6*10^-2])

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',16, 'fontname','Times')
legend('Rotational', 'Divergent','-ve Divergent','Total', 'r^2', 'r^{2/3}', 'location', 'best')
ylabel('<\delta u_i^2> [m^2/s^2]')
xlabel('r [km]')
print('./figures/SF2_GLAD_helm.eps','-depsc', '-r400')

%% SF2 Helm LASER

figure
loglog(GLAD.dist_axis/1e3, LASER.SF2rr, 'color',cols(4,:), 'linewidth',2)
hold all
loglog(GLAD.dist_axis/1e3, LASER.SF2dd, 'color',cols(5,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, -LASER.SF2dd, '--','color',cols(5,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, LASER.SF2rr + LASER.SF2dd, 'color',cols(3,:), 'linewidth',2)

%loglog(GLAD.dist_axis/1e3, GLAD.SF2rr + GLAD.SF2dd, '--', 'color',cols(3,:), 'linewidth',0.8)

ao = 10^-2;
a1 = 10^-2;
xax = [0.001:10:10^3];
yax2 = ao*(xax.^2);
yax1 = 1.3*a1*(xax.^1);
yax23 = 1.4*a1*(xax.^(2/3));
loglog(xax, yax2,'--','color',[0.5 0.5 0.5],'linewidth',0.7)
%loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',0.7)
loglog(xax, yax23,'-','color',[0.5 0.5 0.5],'linewidth',0.7)
%axis([10^-2 10^3 10^-4 6*10^-2])

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',16, 'fontname','Times')
legend('Rotational', 'Divergent','-ve Divergent','Total', 'r^2', 'r^{2/3}', 'location', 'best')
ylabel('<\delta u_i^2> [m^2/s^2]')
xlabel('r [km]')
print('./figures/SF2_LASER_helm_box.eps','-depsc', '-r400')

%% SF2 rotational comparison 

figure
loglog(GLAD.dist_axis/1e3, GLAD.SF2rr, 'color',cols(1,:), 'linewidth',2)
hold all
loglog(GLAD.dist_axis/1e3, LASER.SF2rr, 'color',cols(2,:), 'linewidth',2)

axis([0.01 1000 5e-5 5e-1])

%%
figure
semilogx(GLAD.dist_axis/1e3, LASER.SF2rr./GLAD.SF2rr, 'color',cols(1,:), 'linewidth',2)
hold all 
semilogx(GLAD.dist_axis/1e3, (LASER.SF2rr + LASER.SF2dd)./(GLAD.SF2rr + GLAD.SF2dd),  'color',cols(2,:), 'linewidth',2)

hlines(1, '1b--')

xlabel('r [km]')
ylabel('LASER/GLAD')

legend('Rotational', 'Total', 'location', 'best')

axis([0.01 1000 0 14])
set(gca,'fontsize',16, 'fontname','Times')

print('./figures/SF2_ratios.eps','-depsc', '-r400')


%% SF2 Helmholtz Raw

figure
loglog(GLAD.dist_axis/1e3, GLAD.SF2rr, 'color',cols(1,:), 'linewidth',2)
hold all
loglog(GLAD.dist_axis/1e3, GLAD.SF2dd, 'color',cols(2,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, GLAD.SF2rr + GLAD.SF2dd, 'color',cols(3,:), 'linewidth',2)

loglog(GLAD.dist_axis/1e3, LASER.SF2rr, '--', 'color',cols(1,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, LASER.SF2dd, '--', 'color',cols(2,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, LASER.SF2rr + LASER.SF2dd, '--', 'color',cols(3,:), 'linewidth',2)

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',16, 'fontname','Times')
legend('GLAD rot', 'GLAD div', 'GLAD total', 'LASER rot', 'LASER div', 'LASER total' , 'location', 'best')
ylabel('<\delta u_i^2> [m^2/s^2]')
xlabel('r [km]')
print('./figures/SF2_raw_helm_both.eps','-depsc', '-r400')

%% SF2 Helmholtz Corr

figure
loglog(GLAD.dist_axis/1e3, GLAD.SF2rrcor, 'color',cols(1,:), 'linewidth',2)
hold all
loglog(GLAD.dist_axis/1e3, GLAD.SF2ddcor, 'color',cols(2,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, GLAD.SF2rrcor + GLAD.SF2ddcor, 'color',cols(3,:), 'linewidth',2)

loglog(GLAD.dist_axis/1e3, LASER.SF2rrcor, '--', 'color',cols(1,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, LASER.SF2ddcor, '--', 'color',cols(2,:), 'linewidth',2)
loglog(GLAD.dist_axis/1e3, LASER.SF2rrcor + LASER.SF2ddcor, '--', 'color',cols(3,:), 'linewidth',2)

axis([0.01 1000 5e-5 5e-1])
set(gca,'fontsize',16, 'fontname','Times')
legend('GLAD rot', 'GLAD div', 'GLAD total', 'LASER rot', 'LASER div', 'LASER total' , 'location', 'best')
ylabel('<\delta u_i^2> [m^2/s^2]')
xlabel('r [km]')
print('./figures/SF2_cor_helm_both.eps','-depsc', '-r400')

%% SF3 correction 
GLAD.SF3lllcor = GLAD.SF3lll - 3*GLAD.SF1l.*GLAD.SF2ll + 2*(GLAD.SF1l.^3);
%GLAD.SF3lttcor = GLAD.SF3ltt - 3*GLAD.SF1l.*GLAD.SF2ll + 2*(GLAD.SF1l.^3);

LASER.SF3lllcor = LASER.SF3lll - 3*LASER.SF1l.*LASER.SF2ll + 2*(LASER.SF1l.^3);


%% SF3

figure
loglog(GLAD.dist_axis/1e3, (GLAD.SF3lll + GLAD.SF3ltt), 'color',cols(1,:), 'linewidth',2)
hold all
loglog(LASER.dist_axis/1e3, (LASER.SF3lll + LASER.SF3ltt), 'color',cols(2,:), 'linewidth',2)

axis([0.01 1000  0, 3.5])

a1 = 0.6*10^-4;
xax = [0.001:10:10^3];
yax1 = a1*(xax.^1);
loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',0.7)

loglog(GLAD.dist_axis/1e3, -(GLAD.SF3lll + GLAD.SF3ltt), '--', 'color',cols(1,:), 'linewidth',2)

loglog(LASER.dist_axis/1e3, -(LASER.SF3lll + LASER.SF3ltt), '--','color',cols(2,:), 'linewidth',2)


set(gca,'fontsize',16, 'fontname','Times')
legend('GLAD', 'LASER', 'r^1','location', 'best')
ylabel('$SF3 [m^3/s^3]$', 'interpreter', 'latex')
xlabel('r [km]')

print('./figures/SF3.eps','-depsc', '-r400')


%%
figure
semilogx(GLAD.dist_axis/1e3, (GLAD.SF3lll + GLAD.SF3ltt)./GLAD.dist_axis, 'color',cols(1,:), 'linewidth',2)
hold all
semilogx(LASER.dist_axis/1e3, (LASER.SF3lll + LASER.SF3ltt)./GLAD.dist_axis, 'color',cols(2,:), 'linewidth',2)
grid

%%
figure
semilogx(GLAD.dist_axis/1e3, (GLAD.SF3lll + GLAD.SF3ltt)./GLAD.dist_axis, 'color',cols(1,:), 'linewidth',2)
hold all
semilogx(LASER.dist_axis/1e3, (LASER.SF3lll + LASER.SF3ltt)./GLAD.dist_axis, 'color',cols(2,:), 'linewidth',2)
grid

%%

figure
semilogx(GLAD.dist_axis/1e3, (GLAD.SF3lllcor )./GLAD.dist_axis, 'color',cols(1,:), 'linewidth',2)
hold all
semilogx(GLAD.dist_axis/1e3, (GLAD.SF3lll )./GLAD.dist_axis, '--', 'color',cols(1,:), 'linewidth',2)
semilogx(LASER.dist_axis/1e3, (LASER.SF3lllcor)./GLAD.dist_axis, 'color',cols(2,:), 'linewidth',2)
semilogx(LASER.dist_axis/1e3, (LASER.SF3lll)./GLAD.dist_axis,'--', 'color',cols(2,:), 'linewidth',2)
grid
%% Rotational-divergent decomposition function 
% do the integrals to calculate the decomposition to rotational and divergent part (Lindborg 2015)
% can potentially improve this by using a proper way to calculate the
% integral instead of the nansum
%clear s2rr s2dd

function [SF2rr SF2dd] = helmholtz_decompose(dist_axis, SF2ll, SF2tt)
mid_dist_axis = 0.5*(dist_axis(1:end-1)+dist_axis(2:end));
mid_diff_du = 0.5*((SF2tt(1:end-1) - SF2ll(1:end-1))+(SF2tt(2:end)-SF2ll(2:end)));

    for i =2:length(dist_axis)
        SF2rr(i) = SF2tt(i) + nansum(1./mid_dist_axis(1:i-1).*mid_diff_du(1:i-1).*diff(dist_axis(1:i)));
        SF2dd(i) = SF2ll(i) - nansum(1./mid_dist_axis(1:i-1).*mid_diff_du(1:i-1).*diff(dist_axis(1:i)));
    end
end
