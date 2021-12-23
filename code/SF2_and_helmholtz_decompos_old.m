% This is the code that is used to generate the results of 
% Balwada, LaCasce and Speer 2016
% (not the exact figures, but all the analysis is here).

% Parts of this code are used for the 2021 paper, but are all in separate
% scripts.

% This script was just put here because the codes for 2016 paper were never
% shared anywhere else.

% Code to do the helmholtz decomposition on the observational
% 2nd order structure functions.

%%%%%%
%%%%%% This is old code!!!
%%%%%% Not used in this work!!!
%%%%%%

clear all
close all

traj = load ('../GLAD_data/glad_traj_297.mat');

%%
% calculate time series of pair separation
sep = calculate_seperation_timeseries(traj);

%%
% plevel = [1200 1600];

gamma = 1.5;

dist_bin(1) = 0.1; % in m
dist_bin = gamma.^[0:100]*dist_bin(1);
id = find(dist_bin>800*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));
%     separationr = sqrt(dispersion(d).disp);
%%

clear s2ll s2tt s3lll s3llt s3ltt s3ttt

s2ll  = zeros(length(dist_axis),1);
s2tt  = zeros(length(dist_axis),1);

s3lll = zeros(length(dist_axis),1);
s3ltt = zeros(length(dist_axis),1);
s3llt = zeros(length(dist_axis),1);
s3ttt = zeros(length(dist_axis),1);

npairs = zeros(length(dist_axis),1);
%
% loop for different distance classes
for i =1:length(dist_axis)
    dull = []; dutt = [];
    % loop over different pairs
    for j = 1:length(sep)
        dull_temp = []; dutt_temp = [];
        % find id of the pairs in a particular geographical regime (within a certain distance from each other)
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i));
        
        % loop over the different pairs that lie in the range
        for k =1:length(id)
            % components of vector joining the two particles
            rx(k) = (sep(j).X1(id(k)) - sep(j).X2(id(k)))*cosd(0.5*(sep(j).Y1(id(k))+sep(j).Y2(id(k))));
            ry(k) = (sep(j).Y1(id(k)) - sep(j).Y2(id(k)));
            magr(k) = sqrt(rx(k).^2+ry(k).^2);
            
            % normalize to unit vectors
            rx(k) = rx(k)/magr(k); ry(k) = ry(k)/magr(k);
            
            % components of velocity differences
            dux(k) = (sep(j).U1(id(k))-sep(j).U2(id(k)));
            duy(k) = (sep(j).V1(id(k))-sep(j).V2(id(k)));
            
            % convert to longitudnal and transverse structure functions
            dull_temp(k) = dux(k)*rx(k) + duy(k)*ry(k);
            dutt_temp(k) = duy(k)*rx(k) - dux(k)*ry(k);

        end
        if ~isempty(dull_temp)
            dull = [dull; dull_temp'];
            dutt = [dutt; dutt_temp'];
        end
    end
    
    struct_pairs(i).dull = dull;
    struct_pairs(i).dutt = dutt;
    
    % calculate the second order structure functions
    if ~isempty(dull) & length(dull)>20
        
        s2ll(i) = nanmean(dull.^2);
        s2tt(i) = nanmean(dutt.^2);
        
        s3lll(i) = nanmean(dull.^3);
        s3llt(i) = nanmean(dull.^2.*dutt);
        s3ltt(i) = nanmean(dull.*dutt.^2);
        s3ttt(i) = nanmean(dutt.^3);
        
    else
        s2ll(i) = NaN;
        s2tt(i) = NaN;
        
        s3lll(i) = NaN;
        s3llt(i) = NaN;
        s3ltt(i) = NaN;
        s3ttt(i) = NaN;
    end
    disp(i)
end

%% 

save structure_pairs.mat struct_pairs dist_axis dist_bin -v7.3

%% Calculate the 4th order structure functions
s4l = 0*dist_axis; 
s4t = 0*dist_axis; 
for i =1:length(struct_pairs)
    s4l(i) = nanmean(struct_pairs(i).dull.^4)/(s2ll(i).^2);
    s4t(i) = nanmean(struct_pairs(i).dutt.^4)/(s2tt(i).^2);
end

%% 
return 
%% Plot the 4th order structure functions (measure of intermittency)
figure

loglog(dist_axis/1000, s4l,'linewidth',2), hold all
loglog(dist_axis/1000, s4t,'linewidth',2)
axis([10^-2 10^3 1 20])
ylabel('Flatness')
xlabel('r(km)')
legend('F_l', 'F_t')
grid


%% 3rd order structure functions normalized by r
s3 = s3lll+s3ltt; 

close all 
figure
loglog(dist_axis/1000, abs(s3/dist_axis'),'+-','linewidth', 2), hold all
loglog(dist_axis/1000, -(s3/dist_axis'),'o','linewidth', 2)
loglog(dist_axis/1000, abs(s3lll/dist_axis'),'+-','linewidth', 2), hold all
loglog(dist_axis/1000, -(s3lll/dist_axis'),'o','linewidth', 2)
loglog(dist_axis/1000, abs(s3ltt/dist_axis'),'+-','linewidth', 2), hold all
loglog(dist_axis/1000, -(s3ltt/dist_axis'),'o','linewidth', 2)
legend('|S3_{lll} + S3{ltt}|','-(S3_{lll} + S3{ltt})', '|S3_{lll}|', '-S3_{lll}','|S3{ltt}|','-(S3{ltt})')
grid
axis([10^-2 10^3 10^-13 10^-7])
xlabel('r (km)')
ylabel('(S3)/r')

%% 
colors = get(gca,'ColorOrder');
%% S3t
figure
h(1) = loglog(dist_axis/1000, s3llt+s3ttt,'.-', 'LineWidth',3, 'Color', colors(1,:)) ;
hold all
h(2) = loglog(dist_axis/1000, abs(s3llt),'.-', 'LineWidth',3, 'Color', colors(2,:)) ;
loglog(dist_axis/1000, -(s3llt),'o',  'Color', colors(2,:)) 
h(3) =loglog(dist_axis/1000, abs(s3ttt),'.-', 'LineWidth',3, 'Color', colors(3,:)) ;
loglog(dist_axis/1000, -(s3ttt),'o', 'LineWidth',3, 'Color', colors(3,:)) 


r= loglog(dist_axis/1000, 1e-10*dist_axis.^2);

xlabel('r (km)')
ylabel('S3_t')
set(gca,'FontSize', 20)
axis([0.001 1000 1e-10 1])
legend([h r], {'S3_t', 'S3_{llt}', 'S3_{ttt}','r^2'},'location','northwest')
legend boxoff
saveas(gcf,'S3t_glad.eps','epsc2')
%% S3l
close all
figure
h(1) = loglog(dist_axis/1000, abs(s3lll+s3ltt),'.-', 'LineWidth',3, 'Color', colors(1,:)) ;
hold all
 loglog(dist_axis/1000, -(s3lll+s3ltt),'o', 'Color', colors(1,:)) ;
h(2) = loglog(dist_axis/1000, abs(s3lll),'.-', 'LineWidth',3, 'Color', colors(2,:)) ;
loglog(dist_axis/1000, -(s3lll),'o',  'Color', colors(2,:)) 
h(3) =loglog(dist_axis/1000, abs(s3ltt),'.-', 'LineWidth',3, 'Color', colors(3,:)) ;
loglog(dist_axis/1000, -(s3ltt),'o',  'Color', colors(3,:)) 


r=loglog(dist_axis/1000, 1e-8*dist_axis);

xlabel('r (km)')
ylabel('S3_l')
set(gca,'FontSize', 20)
axis([0.001 1000 1e-10 1])
legend([h r], {'S3_l', 'S3_{lll}', 'S3_{ltt}', 'r'}, 'location','northwest')
legend boxoff
saveas(gcf,'S3l_glad.eps','epsc2')

%% ratios 

figure 
loglog(dist_axis/1000, abs(s3lll)./abs(s3ltt),'.-', 'LineWidth',3)
ylabel('|S3_{lll}|/|S3_{ltt}|')
set(gca,'FontSize', 20)
axis([0.001 1000 0.1 100])
saveas(gcf,'S3l_ratio.eps','epsc2')


%% On a semilog axis
figure
semilogx(dist_axis/1000, s3./dist_axis/2,'o-','linewidth',2)
axis([10^-2 10^3 -1*10^-7 1.5*10^-7])
xlabel('r (km)')
ylabel('<\bf \delta u.\delta u \rm \delta u_l>/2r (m^2/s^3)')
grid on
set(gca,'fontsize',22)

%% 3rd order and 4th order quantities on same plot
figure 

[ax,h1,h2]=plotyy(dist_axis/1000, abs(s3), dist_axis/1000, ...
        s4l,'loglog','loglog')
hold(ax(1))
loglog(ax(1),dist_axis/1000, -s3,'o','linewidth', 2)
hold(ax(2))
loglog(ax(2),dist_axis/1000,s4t,'--','linewidth',2,'color',[ 0.8500    0.3250    0.0980])
set(ax(1),'Ylim',[10^-8 1],'Xlim',[10^-2 10^3],'fontsize',16,'Xaxislocation','top')
set(ax(2),'Ylim',[1 150],'Xlim',[10^-2 10^3],'fontsize',16,'Xaxislocation','top')

a1 = 0.6*10^-4;
xax = [0.001:10:10^3];
yax1 = a1*(xax.^1);
loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',2)

set(h1,'linewidth',2,'marker','+')
set(h2,'linewidth',2)
xlabel('r (km)')
ylabel(ax(1), '<\bf \delta u.\delta u \rm \delta u_l> (m^3/s^3)')
ylabel(ax(2), 'Kurtosis')
legend(ax(2),'longitudinal','transverse')
box on
grid on


%% How many pairs contribute to each separation bin
for i =1:length(struct_pairs)
    npairs(i) = length(struct_pairs(i).dull);
end
%%
figure,
loglog(dist_axis/1000,npairs,'s-','linewidth',2)
%axis([10^-2 10^3 1000 10^8]),
grid
xlabel('r (km)')
ylabel('Number of pairs')

%% do the integrals to calculate the decomposition to rotational and divergent part (Lindborg 2015)
% can potentially improve this by using a proper way to calculate the
% integral instead of the nansum
clear s2rr s2dd
mid_dist_axis = 0.5*(dist_axis(1:end-1)+dist_axis(2:end));
mid_diff_du = 0.5*((s2tt(1:end-1)-s2ll(1:end-1))+(s2tt(2:end)-s2ll(2:end)));

for i =2:length(dist_axis)
    s2rr(i) = s2tt(i) + nansum(1./mid_dist_axis(1:i-1).*mid_diff_du(1:i-1)'.*diff(dist_axis(1:i)));
    s2dd(i) = s2ll(i) - nansum(1./mid_dist_axis(1:i-1).*mid_diff_du(1:i-1)'.*diff(dist_axis(1:i)));
end

%% plot the components of 2nd order structure functions
% close all
figure
loglog(dist_axis/1000, s2ll,'linewidth',2)
hold all
loglog(dist_axis/1000, s2tt,'linewidth',2)
hold all
loglog(dist_axis/1000, s2tt+s2ll, '-','linewidth',2)
hold all
% loglog(dist_axis/1000, s2rr+s2dd, '+','linewidth',2)
loglog(dist_axis/1000, s2rr, '-','linewidth',2)
loglog(dist_axis/1000, s2dd, '-','linewidth',2)

ao = 10^-2;
a1 = 10^-2;
xax = [0.001:10:10^3];
yax2 = ao*(xax.^2);
yax1 = 1.3*a1*(xax.^1);
yax23 = 1.4*a1*(xax.^(2/3));
loglog(xax, yax2,'--','color',[0.5 0.5 0.5],'linewidth',2)
loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',2)
loglog(xax, yax23,'-','color',[0.5 0.5 0.5],'linewidth',2)
axis([10^-2 10^3 10^-4 6*10^-2])
xlabel('r (km)')
ylabel('<\delta u^2> (m^2/s^2)')
h=legend('D_{ll}','D_{tt}','D_{tt} + D_{ll}','D_{rr}','D_{dd}')
set(h,'fontsize',16)
set(gca,'fontsize',16)
grid


