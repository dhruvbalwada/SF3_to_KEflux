%% Plot of trajectories for paper 
% Dhruv Balwada, 2 Feb 2021

clear all
close all

traj = load ('/Users/dhruvbalwada/Google Drive/Main_projects/GOM/GLAD_data/glad_traj_297.mat');

%% 
close all

ndays = 10; 
figure
for i = 1:size(traj.X, 2) 
    id_start(i) = find(~isnan(traj.X(:,i)), 1);
    
    plot(traj.X(id_start(i), i), traj.Y(id_start(i), i), 'k.', 'markersize',12)
    plot(traj.X(id_start(i):id_start(i)+4*24*ndays, i), traj.Y(id_start(i):id_start(i)+4*24*ndays, i))
    hold all
end

axis([-90 -86 26, 30])
xlabel('Longitude','interpreter','latex')
ylabel('Latitude','interpreter','latex')
set(gca,'FontSize',18,'FontName','Times')


print('traj_plot.eps','-depsc', '-r400')