

%% 
% get the color axis
colors = get(gca,'ColorOrder');
close all


%% %% Energy injections  %% %%
% There are two ways to show the energy injection, either in regular form
% or in variance preserving form. The variance preserving form changes the
% height of the curves to account for the non-linear squeezing that takes
% place on a log axis. This the more appropriate form visually representing
% quantities for which the area under the curve is important. 

close all
%% Plot of energy injection on k axis (var preserving)
figure
%shadedErrorBar_semilogx((kf)*1e3, kf.*mean_ebs(1:end-1), CI_ebs(:,1:end-1) ...
%                ,{'o-','linewidth',2,'color', colors(1,:)}, 0.6)
%hold all 
semilogx((kf)*1e3, ebs(1:end-1,:).*kf', 'linewidth',.5)%,'color', colors(2,:))
hold all
semilogx((kf)*1e3, mean_ebs(1:end-1).*kf', 'linewidth',2,'color', 'k')

semilogx(1e-3, mean_ebs(end), '.', 'color', 'b', 'Markersize', 30)

%semilogx((1./kf)/1e3, mean_ebs(1:end-1).*(kf)', 'linewidth',2,'color', colors(2,:))
%axis([10^-2 10^3 0 10^-7])

%legend('Errorbars: 5th and 95th percentiles', 'Mean', 'location','northeast')
grid on
xlabel('$$k [km^{-1}]$$','interpreter','latex')
ylabel('$$k*\epsilon_j [m^2s^{-3}]$$','interpreter','latex')
%title('Energy Input','interpreter','latex')
set(gca,'FontSize',18,'FontName','Times')

print(['energy_inject_', experiment,'_', inv_style, '.png'],'-dpng', '-r400')
%
%if strcmp(experiment, 'LASER')
%    print(['energy_inject_', experiment, inv_style, '.png'],'-dpng', '-r400')
%else
%    print('energy_inject_GLAD_LS.png','-dpng', '-r400')
%end

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

%% Plot of 3rd Order SF
figure 
loglog(r/1e3, abs(mean_S),'-' , 'color',colors(1,:),'linewidth',2,'markersize',7)
hold on
loglog(r/1e3,mean_S,'+','color',colors(1,:),'linewidth',2,'markersize',7)
loglog(r/1e3,-mean_S,'o','color',colors(2,:),'linewidth',2,'markersize',7)

loglog(R/1e3,mean_Vt,'r','linewidth',2,'color',colors(5,:))
loglog(R/1e3,-mean_Vt,'r--','linewidth',2, 'color',colors(5,:))

axis([10^-2 10^3 10^-7 1])
xticks([10^-2, 10^-1, 1, 10, 100, 1000])

leg=legend('$|V|$', '$+ V$','$- V$','$$+V_{fit}$$','$$-V_{fit}$$');
set(leg,'interpreter','latex', 'location','best')
xlabel('$$r [km]$$','interpreter','latex')
ylabel('$$V (r) [m^3s^{-3}]$$','interpreter','latex')
set(gca,'FontSize',20,'FontName','Times')

print(['S3_and_fit_', experiment,'_', inv_style, '.png'],'-dpng', '-r400')

%
%if strcmp(experiment, 'LASER')
%    print('S3_and_fit_LASER_LS.png','-dpng', '-r400')
%else
%    print('S3_and_fit_GLAD_LS.png','-dpng', '-r400')
%end

%% Spectral Flux plot 
figure
shadedErrorBar_semilogx(kf*1e3, mean_SpecFlux, CI_SpecFlux ...
                ,{'o-','linewidth',2,'color', colors(1,:)}, 0.6)
hold all 
%semilogx(kf*1e3, mean_SpecFlux, 'linewidth',2,'color', colors(2,:))
%grid on
%leg=legend('Errorbars: 5th and 95th percentiles','Median','Mean');
%set(leg,'interpreter','latex', 'location','best')
yline(0, '--')
xlabel('$$k [1/km]$$','interpreter','latex')
ylabel('$$F(k) [m^2 s^{-3}]$$','interpreter','latex')
set(gca,'FontSize',20,'FontName','Times')

%axis([10^-3 100 -18*10^-8 10*10^-8])
xticks([10^-3, 10^-2, .1, 1, 10, 100])

print(['spectral_flux_', experiment,'_', inv_style, '.png'],'-dpng', '-r400')

%if strcmp(experiment, 'LASER')
%    print('spectral_flux_LASER_LS.png','-dpng', '-r400')
%else
%    print('spectral_flux_GLAD_LS.png','-dpng', '-r400')
%end