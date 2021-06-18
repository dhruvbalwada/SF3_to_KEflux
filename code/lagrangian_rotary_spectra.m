% Dhruv Balwada 
% 6 May 2021

% code to compute the lagrangian spectral

clear all
close all 

exp = 'LASER'; 

if strcmp(exp, 'GLAD')
    load ~/work_root/GoMexico_drifters/GLAD_15min_filtered/traj_struct_GLAD_15min_03_May_2021.mat
    load ../data/traj_depth_GLAD.mat
else
    load ~/work_root/GoMexico_drifters/LASER_SPOT_15min_filtered/traj_structs_LASER_15min_03_May_2021.mat
    load ../data/traj_depth_LASER.mat
    drifter = drifterL;
end
%%

% npts in 30days with 15min resolution
n30 = 28*24*4; 

%%

npts = zeros(length(drifter), 1); 

for i = 1:length(drifter)
    npts(i) = length(drifter(i).time); 
end

%% See how many trajectories have more than
figure
hist(npts, 50)
vlines([n30, 2*n30, 3*n30])

%% Break up into little bits of uniform size
% This strategy leaves out data at the end of the trajectory. 
% But maybe it is practically not very bad. 

nseg = 1; 
clear CV

for i = 1:length(drifter) 
    nparts = npts(i)/n30;
    for j = 1:floor(nparts)
        CV(1:n30, nseg) = drifter(i).u((j-1)*n30+1: j*n30) + sqrt(-1)*drifter(i).v((j-1)*n30+1: j*n30);
        %hseg = 
        nseg = nseg + 1;
    end
end

%%
psi = sleptap(n30); 
%%
%[F, SPP, SNN, SPN] = mspec(15*60, CV, psi, 'cyclic');
[F, SPP, SNN, SPN] = mspec(15*60, CV, [], 'cyclic');

%
f_coriollis = 2 * (1/24/3600) * sind(28); % this does not have a 2pi as it is in units of cycles/day

%%
figure
loglog(F/f_coriollis, mean(SPP,2) )
hold all 
loglog(F/f_coriollis, mean(SNN,2) )

loglog(F/f_coriollis, 1e-8*F.^(-2), '--', 'color','k' )

legend('Cyclonic', 'Anticyclonic', '\omega^{-2}')
axis([5e-2 10 1 1e5])

xlabel('\omega/f')
ylabel('$\hat{E} (\omega)$', 'interpreter','latex')
set(gca,'FontSize',20,'FontName','Times')

%%
if strcmp(exp,'LASER')
    print('./figures/lagr_spec_LASER.eps','-depsc', '-r400')
else
    print('./figures/lagr_spec_GLAD.eps','-depsc', '-r400')
end
