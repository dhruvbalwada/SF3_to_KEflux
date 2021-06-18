%% 3rd Order Structure Function to Kinetic Energy Spectral Flux
% Dhruv Balwada, 29 January 2021
% 
% 

clear all 
close all 


%%

experiment = 'GLAD';

if strcmp(experiment, 'LASER')
    load ../data/LASER_S3_deep500_boot_strap.mat
else
    load ../data/GLAD_S3_deep500_boot_strap.mat
end
%%

nsamps = size(s3lll,2);

r=dist_axis;
Nr=length(r);

% pick data points larger than 10m for fitting
ns=find(dist_axis>=20 ,1);
ne=find(dist_axis<=500e3 ,1,'last');

R=r(ns:ne);

NR=length(R);
lf=R;

%%

%%
% kf=10.^(-5:0.1:-1);
% kf=1e-3:5e-3:1e-1;

kf=1./lf;
Nk=length(kf);

%e=zeros(1,Nk+1);
ebs = zeros(NR+1,nsamps);
S = zeros(Nr, nsamps);
Vt = zeros(NR,nsamps);
SpecFlux = zeros(Nk, nsamps); 

for n=1:nsamps
    
    S(:,n) =s3lll(:,n)' +s3ltt(:,n)';

    V=S(ns:ne, n)';
    
    % build the matrix A
    A=zeros(NR,Nk+1);
    for j=1:Nk
        A(:,j)=-4/kf(j)*besselj(1,kf(j)*R)';
    end
    A(:,end)=2*R';
    
    % normalize the magnitude
    
    for j=1:NR
        A(j,:)=A(j,:)./abs(R(j));
    end
    V=V./abs(R);
    
    % estimate the epsilons using least squares
    ebs(:,n) = lsqnonneg(A,V');
    
    % What does the V look like
    Vt(:,n) =2*ebs(end,n)*R;
    for j=1:Nk
        Vt(:,n)=Vt(:,n)-4*ebs(j,n)./kf(j).*besselj(1,kf(j).*R)';
    end
    
    % Convert V to spectral fluxes 
    SpecFlux(1,n) = -ebs(end,n); 
    for j = 2:Nk
        SpecFlux(j,n) = SpecFlux(j-1,n) + ebs(end-j + 1,n);
    end
    
end

SpecFlux = flipud(SpecFlux); 
%% Estimate confidence intervals 

clear CI_ebs CI_Vt

for i = 1:size(ebs,1)
    CI_ebs(:,i) = prctile(ebs(i,:), [95, 5]); 
end
for i = 1:size(Vt,1)
    CI_Vt(:,i) = prctile(Vt(i,:), [95,5]);
    CI_SpecFlux(:,i) = prctile(SpecFlux(i,:), [95,5]);
end

median_ebs = median(ebs,2); 
mean_ebs = nanmean(ebs,2); 
median_Vt = median(Vt,2);
mean_Vt = nanmean(Vt,2);

median_S = median(S,2);
mean_S = nanmean(S,2);

median_SpecFlux = median(SpecFlux,2);
mean_SpecFlux = nanmean(SpecFlux,2);


%% 
% get the color axis
colors = get(gca,'ColorOrder');

%% Plot of energy injection 
%close all 
figure
shadedErrorBar_semilogx(R/1e3, mean_ebs(1:end-1), CI_ebs(:,1:end-1) ...
                ,{'o-','linewidth',2,'color', colors(1,:)}, 0.6)
%hold all 
%semilogx(R/1e3, mean_ebs(1:end-1), 'linewidth',2,'color', colors(2,:))
axis([10^-2 10^3 0 10^-7])

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

%% Plot of energy injection 
%close all 
figure
neg = mean_ebs - CI_ebs(2,:)';
pos = - mean_ebs + CI_ebs(1,:)';
errorbar(R/1e3, mean_ebs(1:end-1), neg(1:end-1), pos(1:end-1), '*', 'linewidth',2); 
hold all 
set(gca,'XScale','log')
axis([10^-2 10^3 0 10^-7])
xticks([10^-2, 10^-1, 1, 10, 100, 1000])

%legend('Errorbars: 5th and 95th percentiles', 'Median', 'Mean', 'location','northeast')

xlabel('$$r [km]$$','interpreter','latex')
ylabel('$$\epsilon_j [m^2s^{-3}]$$','interpreter','latex')
%title('Energy Input','interpreter','latex')
set(gca,'FontSize',20,'FontName','Times')

%print('energy_input_2.eps','-depsc', '-r400')


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

%
if strcmp(experiment, 'LASER')

    print('S3_and_fit_LASER.png','-dpng', '-r400')
else
    print('S3_and_fit_GLAD.png','-dpng', '-r400')
end

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

axis([10^-3 100 -18*10^-8 10*10^-8])
xticks([10^-3, 10^-2, .1, 1, 10, 100])

if strcmp(experiment, 'LASER')

    print('spectral_flux_LASER.png','-dpng', '-r400')
else
    print('spectral_flux_GLAD.png','-dpng', '-r400')
end