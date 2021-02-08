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


%print('traj_plot.eps','-depsc', '-r400')

%% 
close all

ndays = 10;
figure
m_proj('Lambert','lon',[-91 -84.5],'lat',[26 31]); 

%[elev, long, latg] = m_elev([-93 -82 26 31]);
[elev, latg, long] = mygrid_sand2([-91.5 -83.5 26 31]);

clear CS CH
[CS, CH] = m_contourf(long-360, latg, elev,  'edgecolor','none');

hold all
%[CS,CH]=m_elev('contourf',[-3500:250:0 ],'edgecolor','none');

m_grid('linestyle','none','tickdir','in','linewidth',1,'FontSize',20,'FontName','Times','interpreter','latex');

%
col = m_colmap('blues',80);
colormap([ col;[0 0 0] ]);
brighten(.5);


%m_coast('speckle',[0 .6 0])
%m_coast('patch',[0 .0 0])
% m_gshhs_i('color','k');
ax=m_contfbar(1.02,[.5 .8],CS,CH);

%
title(ax,{'Depth [m]',''}); % Move up by inserting a blank line

for i = 1:size(traj.X, 2) 
    id_start(i) = find(~isnan(traj.X(:,i)), 1);
    
    m_plot(traj.X(id_start(i), i), traj.Y(id_start(i), i), 'k.', 'markersize',12)
    m_plot(traj.X(id_start(i):id_start(i)+4*24*ndays, i), traj.Y(id_start(i):id_start(i)+4*24*ndays, i), 'linewidth',1.5)
    hold all
end

print('traj_plot.eps','-depsc', '-r400')
%set(gca,'FontSize',20,'FontName','Times')