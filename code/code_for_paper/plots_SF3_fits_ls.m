%% Plots of SF3 and related metrics
% Dhruv Balwada
% 18 December 2021

%% 
% get the color axis
colors = get(gca,'ColorOrder');
close all

if strcmp(experiment, 'GLAD')
    col_num = 7;
else
    col_num = 1;
end
%% %% Energy injections  %% %%
% There are two ways to show the energy injection, either in regular form
% or in variance preserving form. The variance preserving form changes the
% height of the curves to account for the non-linear squeezing that takes
% place on a log axis. This the more appropriate form visually representing
% quantities for which the area under the curve is important. 

%% Plot of energy injection on k axis (var preserving)
figure
shadedErrorBar_semilogx((kf)*1e3, mean_ebs(1:end-1), CI_ebs(:,1:end-1) ...
                ,{'o-','linewidth',2,'color', colors(col_num,:),'Markersize',5}, 0.6)

%grid on
xlabel('$$k [1/km]$$','interpreter','latex')
ylabel('$$k.\epsilon_j [m^2s^{-3}]$$','interpreter','latex')
%title('Energy Input','interpreter','latex')
set(gca,'FontSize',18,'FontName','Times')
axis([10^-3 100 -2*10^-6 2*10^-6])

xticks([10^-3, 10^-2, .1, 1, 10, 100])

ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'xscale','log');
set(ax2, 'XTick', get(ax1, 'XTick') , 'YTick', get(ax1, 'YTick'));
%OppTickLabels = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k'}
OppTickLabels = {'' '100km' '10km' '1km' '100m' '10m' }
% Set the x-tick and y-tick  labels for the second axes
set(ax2, 'XTickLabel', OppTickLabels,'YTickLabel',{}, 'FontSize',18,'FontName','Times')


print(['energy_inject_', experiment,'_', inv_style, '.png'],'-dpng', '-r400')


%% Spectral Flux plot 
figure
shadedErrorBar_semilogx(kf*1e3, mean_SpecFlux, CI_SpecFlux ...
                ,{'o-','linewidth',2,'color', colors(col_num,:),'Markersize',5}, 0.6)
hold all 
%semilogx(kf*1e3, mean_SpecFlux, 'linewidth',2,'color', colors(2,:))
%grid on
%leg=legend('Errorbars: 5th and 95th percentiles','Median','Mean');
%set(leg,'interpreter','latex', 'location','best')
yline(0, '--')
xlabel('$$k [1/km]$$','interpreter','latex')
ylabel('$$F(k) [m^2 s^{-3}]$$','interpreter','latex')
set(gca,'FontSize',20,'FontName','Times')

axis([10^-3 100 -60*10^-8 60*10^-8])
xticks([10^-3, 10^-2, .1, 1, 10, 100])


ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'xscale','log');
set(ax2, 'XTick', get(ax1, 'XTick') , 'YTick', get(ax1, 'YTick'));
%OppTickLabels = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k'}
OppTickLabels = {'' '100km' '10km' '1km' '100m' '10m' }
% Set the x-tick and y-tick  labels for the second axes
set(ax2, 'XTickLabel', OppTickLabels,'YTickLabel',{}, 'FontSize',18,'FontName','Times')
print(['spectral_flux_', experiment,'_', inv_style, '.png'],'-dpng', '-r400')


%% SF3 with errorbars 
CI_SF3_by_r = CI_SF3'./r';
%

figure
shadedErrorBar_semilogx(r/1e3, mean_SF3./r', CI_SF3_by_r ...
                ,{'o-','linewidth',2,'color', colors(col_num,:),'Markersize',5}, 0.6)

hold all
plot(R/1e3, (mean_Vt)./R','--', 'linewidth',2, 'color','k')

yline(0, '--')


%grid on
xlabel('$$r [km]$$','interpreter','latex')
ylabel('$$ V(r)/r$$   $$[m^2s^{-3}]$$','interpreter','latex')
%title('Energy Input','interpreter','latex')
set(gca,'FontSize',18,'FontName','Times')
axis([10^-2 500 -3*10^-7 6*10^-7])

print(['Vbyr_', experiment,'_', inv_style, '.png'],'-dpng', '-r400')

%%
stop

%%
LASER_SF3 = load('../data/LASER_S3_deep500_box_constrained_block_boot_strap_Ldof.mat');
GLAD_SF3 = load('../data/GLAD_S3_deep500_block_boot_strap_Ldof.mat');

GLAD_SF3 = nanmean(GLAD_SF3.SF3, 2); 
LASER_SF3 = nanmean(LASER_SF3.SF3, 2); 

%%

figure
G = loglog(r/1e3, GLAD_SF3,'-','color',colors(7,:),'linewidth',2,'markersize',7);
hold on
loglog(r/1e3,-GLAD_SF3,'--','color',colors(7,:),'linewidth',2,'markersize',7)

L = loglog(r/1e3, LASER_SF3,'-','color',colors(1,:),'linewidth',2,'markersize',7);
loglog(r/1e3,-LASER_SF3,'--','color',colors(1,:),'linewidth',2,'markersize',7)

lin = loglog(r/1e3, 1e-7*r, '--', 'color', [0.5, 0.5, 0.5]);

axis([10^-2 10^3 10^-7 0.3])
xticks([10^-2, 10^-1, 1, 10, 100, 1000])

legend([G,L, lin], 'GLAD', 'LASER','$r^1$', 'location', 'best','interpreter','latex')
set(gca,'FontSize',22,'FontName','Times')
xlabel('$$r [km]$$','interpreter','latex')
ylabel('$$V (r) [m^3s^{-3}]$$','interpreter','latex')

print(['S3_both.png'],'-dpng', '-r400')
%print(['S3_both.eps'],'-deps2')

%% Plot of 3rd Order SF and fit
figure 
loglog(r/1e3, abs(mean_SF3),'-' , 'color',colors(1,:),'linewidth',2,'markersize',7)
hold on
loglog(r/1e3,mean_SF3,'+','color',colors(1,:),'linewidth',2,'markersize',7)
loglog(r/1e3,-mean_SF3,'o','color',colors(2,:),'linewidth',2,'markersize',7)

loglog(R/1e3,mean_Vt,'r','linewidth',2,'color',colors(5,:))
loglog(R/1e3,-mean_Vt,'r--','linewidth',2, 'color',colors(5,:))

axis([10^-2 10^3 10^-7 1])
xticks([10^-2, 10^-1, 1, 10, 100, 1000])

leg=legend('$|V|$', '$+ V$','$- V$','$$+V_{fit}$$','$$-V_{fit}$$');
set(leg,'interpreter','latex', 'location','best')
xlabel('$$r [km]$$','interpreter','latex')
ylabel('$$V (r) [m^3s^{-3}]$$','interpreter','latex')
set(gca,'FontSize',20,'FontName','Times')

%print(['S3_and_fit_', experiment,'_', inv_style, '.png'],'-dpng', '-r400')

%% SF3 with errorbars 

figure 
plot(r/1e3, (mean_SF3)./r')
hold all
plot(R/1e3, (mean_Vt)./R')
plot(r/1e3, (CI_SF3(1,:)')./r')
plot(r/1e3, (CI_SF3(2,:)')./r')

set(gca, 'XScale','log', 'YScale','linear')

axis([10^-2 1000 -5*10^-7 5*10^-7])

grid


%% 
% Attempt to flip and stuff the errorbar, but baat nahi bani. Zero ke paas
% kaam nahi kar rha hai. 

CI_SF3_flipd = CI_SF3; 
CI_SF3_flipd(:,mean_SF3>0) =  flipud(CI_SF3(:,mean_SF3>0) );

%%

figure 
plot(r/1e3, abs(mean_SF3))
hold all
plot(r/1e3, abs(CI_SF3_flipd(1,:)'))
plot(r/1e3, abs(CI_SF3_flipd(2,:)'))

set(gca, 'XScale','log', 'YScale','log')
grid

%%%%%%%%%%%%%%
%% Other plots %%
%%%%%%%%%%%%%%

stop 
%% Plot of energy injection on r axis (var preserving)

plot_flag = 0 ;

if plot_flag==1
figure
%shadedErrorBar_semilogx((kf)*1e3, kf.*mean_ebs(1:end-1), CI_ebs(:,1:end-1) ...
%                ,{'o-','linewidth',2,'color', colors(1,:)}, 0.6)
%hold all 
semilogx((1./kf)/1e3, mean_ebs(1:end-1).*kf', 'linewidth',2,'color', colors(2,:))

%semilogx((1./kf)/1e3, mean_ebs(1:end-1).*(kf)', 'linewidth',2,'color', colors(2,:))
%axis([10^-2 10^3 0 10^-7])

%legend('Errorbars: 5th and 95th percentiles', 'Mean', 'location','northeast')

xlabel('$$k [km^{-1}]$$','interpreter','latex')
ylabel('$$k*\epsilon_j [m^2s^{-3}]$$','interpreter','latex')
%title('Energy Input','interpreter','latex')
set(gca,'FontSize',18,'FontName','Times')

% Plot of energy injection on r axis (not var preserving)
figure
shadedErrorBar_semilogx((1./kf)/1e3, mean_ebs(1:end-1), CI_ebs(:,1:end-1) ...
                ,{'o-','linewidth',2,'color', colors(1,:)}, 0.6)
%hold all 
%semilogx((1./kf)/1e3, mean_ebs(1:end-1)./dk', 'linewidth',2,'color', colors(2,:))
%axis([10^-2 10^3 0 10^-7])

%legend('Errorbars: 5th and 95th percentiles', 'Mean', 'location','northeast')

xlabel('$$r [km]$$','interpreter','latex')
ylabel('$$\epsilon_j [m^2s^{-3}]$$','interpreter','latex')
%title('Energy Input','interpreter','latex')
set(gca,'FontSize',18,'FontName','Times')
%
if strcmp(experiment, 'LASER')

    print('energy_inject_LASER.png','-dpng', '-r400')
else
    print('energy_inject_GLAD.png','-dpng', '-r400')
end
%print('energy_input.eps','-depsc', '-r400')
end
%% Plot of energy injection 
plot_errbar_flag = 0 ;

if plot_errbar_flag==1
%close all 
figure
neg = mean_ebs - CI_ebs(2,:)';
pos = - mean_ebs + CI_ebs(1,:)';
errorbar((1./kf)/1e3, mean_ebs(1:end-1).*kf', neg(1:end-1).*kf', pos(1:end-1).*kf', '*', 'linewidth',2); 
hold all 
set(gca,'XScale','log')
%axis([10^-2 10^3 0 10^-7])
xticks([10^-2, 10^-1, 1, 10, 100, 1000])

%legend('Errorbars: 5th and 95th percentiles', 'Median', 'Mean', 'location','northeast')

xlabel('$$r [km]$$','interpreter','latex')
ylabel('$$\epsilon_j [m^2s^{-3}]$$','interpreter','latex')
%title('Energy Input','interpreter','latex')
set(gca,'FontSize',20,'FontName','Times')

%print('energy_input_2.eps','-depsc', '-r400')
end
