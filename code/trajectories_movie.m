%% Movies of trajectories
% Dhruv Balwada, 13 May 2021

clear all
close all

traj = load ('/Users/dhruvbalwada/work_root/GoMexico_drifters/GLAD_15min_filtered/traj_mat_GLAD_15min_04_May_2021.mat');
%traj = load ('/Users/dhruvbalwada/work_root/GoMexico_drifters/LASER_SPOT_15min_filtered/traj_mat_LASER_15min_04_May_2021.mat');

%%
traj.X = traj.trajmat_X;
traj.Y = traj.trajmat_Y;
traj.U = traj.trajmat_U;
traj.V = traj.trajmat_V;

%% Get bathy
[elev, latg, long] = mygrid_sand2([-99 -79 17 33]);
%%
%close all

figure
% setup projection
m_proj('Lambert','lon',[-93 -83],'lat',[22.5 31]);

% plot bathy
clear CS CH
[CS, CH] = m_contourf(long-360, latg, elev, [-7500:250:0 ], 'edgecolor','none');
hold all

[CS1, CH1] = m_contour(long-360, latg, elev, [-500:-499], '--', 'edgecolor','k');

m_grid('linestyle','none','tickdir','in','linewidth',1,'FontSize',20,'FontName','Times','interpreter','latex');

%
col = m_colmap('blues',80);
%colormap([ col ]);
colormap([ col;[0 0 0] ]);
brighten(.5);
caxis([-4000 0])
%

%m_coast('speckle',[0 .6 0])
%m_coast('patch',[0 .0 0])
% m_gshhs_i('color','k');
ax=m_contfbar(1.02,[.5 .8],CS,CH);

%
title(ax,{'Depth [m]',''}); % Move up by inserting a blank line

%print('temp.eps','-depsc', '-r400')
%print('temp.png','-dpng', '-r400')
%
ndays=20;
nstep = ndays*4*24;
for i = 1:size(traj.X, 2)
    id_start(i) = 1; %find(~isnan(traj.X(:,i)), 1);
    
    %m_plot(traj.X(id_start(i), i), traj.Y(id_start(i), i), 'k.', 'markersize',12)
    %m_plot(traj.X(id_start(i):id_start(i)+4*24*ndays, i), traj.Y(id_start(i):id_start(i)+4*24*ndays, i), 'linewidth',1.1)
    m_plot(traj.X(id_start(i)+4*24*ndays:id_start(i)+4*24*(ndays+10), i), traj.Y(id_start(i)+4*24*ndays:id_start(i)+4*24*(ndays+10), i), 'linewidth',1.1)
    hold all
end

%
for i = 1:size(traj.X, 2)
    id_start(i) = find(~isnan(traj.X(:,i)), 1);
    
    m_plot(traj.X(id_start(i), i), traj.Y(id_start(i), i), 'k.', 'markersize',7)
    m_plot(traj.X(1+4*24*ndays, i), traj.Y(1+4*24*ndays, i), 'r.', 'markersize',7)
    %m_plot(traj.X(id_start(i):id_start(i)+4*24*ndays, i), traj.Y(id_start(i):id_start(i)+4*24*ndays, i), 'linewidth',1.5)
    %hold all
end

title(datestr(traj.T_axis(nstep)))

print('temp1.png','-dpng', '-r400')

%% Find start and end ids for each drifter

for i = 1:size(traj.X, 2)
    id_start(i) = find(~isnan(traj.X(:,i)), 1);
    id_end(i) = find(~isnan(traj.X(:,i)), 1, 'last');
end

%% Setup some parameters
close all 

dtres = 0.25; % resolution of data (hours)
nres = floor(1/dtres); % num of frames per hour

dt = 3*nres; % [hours*nframes] the number is the number of hours we want to plot in each frame

Ndt = 5 * 24 *nres; % [days * hours_perday * num frames per hour] the number is the number of days

%%

for i = 1:floor(length(traj.T_axis)/dt) % loop over time (will make frames)
    close all
    figure
    m_proj('Lambert','lon',[-93 -83],'lat',[22.5 31]);
    
    % plot bathy
    clear CS CH
    [CS, CH] = m_contourf(long-360, latg, elev, [-7500:250:0 ], 'edgecolor','none');
    hold all
    [CS1, CH1] = m_contour(long-360, latg, elev, [-500:-499], '--', 'edgecolor','k');
    m_grid('linestyle','none','tickdir','in','linewidth',1,'FontSize',20,'FontName','Times','interpreter','latex');
    col = m_colmap('blues',80);
    colormap([ col;[0. 0. 0.] ]);
    brighten(.5);
    caxis([-4000 0])
    ax=m_contfbar( 0.98,[.5 .8],CS,CH, 'FontName','Times', 'FontSize', 14);
    
    % plot trajectory
    
    % figure out what times to plot
    idplot2 = i*dt;
    idplot1 = 1; % keep the first plotting point at the start if we are initial phase
    if idplot2 > Ndt
       idplot1 = idplot2 - Ndt +1 ;  % drag the first plotting point along if we are plotting long
    end
    
    % plot the track 
    for k = 1:size(traj.X, 2) % loop over trajectories
        m_plot(traj.X(idplot1:idplot2, k), traj.Y(idplot1:idplot2, k), 'linewidth',1.1)
    end
    
    % plot the start and ed points
    num=0;
    for k = 1:size(traj.X, 2) % loop over trajectories
        if (idplot1<= id_start(k) & idplot2>= id_start(k)) % starting point is b/w the range being plotted
            m_plot(traj.X(id_start(k), k), traj.Y(id_start(k), k), 'k.', 'markersize',7) % mark start point as black
        %elseif idplot1> id_start(k) % starting point has been left behind
            %m_plot(traj.X(idplot1, k), traj.Y(idplot1, k), 'g.', 'markersize',7) % mark this initial moving point as green
        end
        
        
        if (idplot2<= id_end(k) & idplot2>= id_start(k)) 
            m_plot(traj.X(idplot2, k), traj.Y(idplot2, k), 'r.', 'markersize',7) % mark this end moving point as red
            num=num+1; % count number of moving drifters
        elseif (idplot2>= id_start(k) & idplot1<id_end(k))
            m_plot(traj.X(id_end(k), k), traj.Y(id_end(k), k), 'c.', 'markersize',7) 
        end
        
    end

    %m_text(-92, 30.5, num2str(num), 'color',[1, 1, 1], 'Fontsize', 16,'Fontweight', 'bold')
    
    title([datestr(traj.T_axis(idplot2)), ', ', num2str(num), ' drifters'],'Fontsize', 16, 'FontName','Times')
    
    %print(['./movie_LASER/frame', num2str(i, '%05.f'),'.png'],'-dpng', '-r400')
    print(['./movie_GLAD/frame', num2str(i, '%05.f'),'.png'],'-dpng', '-r400')

end

