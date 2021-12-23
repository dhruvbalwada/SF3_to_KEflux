%% Test to see if vectorized code works
% Dhruv Balwada
% May 2021

clear all 

%% Make fake data

X = [1,2,3]';
Y = [10, 12, 15]';

U = [0.1, 0.2, 0.3]';
V = [-0.1, -0.4, -0.3]';

%X = 10*rand(100,1);
%Y = 10*rand(100,1);
%U = rand(100,1);
%V = rand(100,1);

%% Convert to traj structure form
traj.X = X';
traj.Y = Y';
traj.U = U';
traj.V = V';

%% Define separation bins

gamma = 1.5;

dist_bin(1) = 0.1; % in m
dist_bin = gamma.^[0:100]*dist_bin(1);

% get rid of bins that are too large
id = find(dist_bin>800*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

%% Old Way
tic
sep = calculate_seperation_timeseries(traj);

% Bin the separation data 

% loop for different distance classes
% this is quite slow.
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
    
    %disp(i)
end

toc

%% New way

n = 1;
for i = 1:2
    for j=i+1:3
        rx_test(n) = (X(i) - X(j))*cosd(0.5* (Y(i)+Y(j)));
        ry_test(n) = (Y(i) - Y(j));
        magr_test(n) = sqrt(rx_test(n)^2 + ry_test(n)^2);
        
        rx_test(n) = rx_test(n)/magr_test(n); 
        ry_test(n) = ry_test(n)/magr_test(n); 
        
        dux_test(n) = U(i) - U(j);
        duy_test(n) = V(i) - V(j);
        
        dul_test(n) = dux_test(n)*rx_test(n) + duy_test(n)*ry_test(n); 
        dut_test(n) = duy_test(n)*rx_test(n) - dux_test(n)*ry_test(n); 
        
        n=n+1;
    end
end

%%
% Input data, specified as a numeric matrix of size m-by-n. Rows correspond
% to individual observations, and columns correspond to individual variables.

tic
Xvec = [X, Y];

dist = pdist(Xvec, @dist_geo);

rx_new = pdist(Xvec, @dist_rx); 
ry_new = pdist(Xvec, @dist_ry); 

magr_new = sqrt(rx_new.^2 + ry_new.^2);

rx_new = rx_new./magr_new; ry_new = ry_new./magr_new;

dux_new = pdist(U, @dist_du);
duy_new = pdist(V, @dist_du);

dull_new = dux_new.*rx_new + duy_new.*ry_new;
dutt_new = duy_new.*rx_new - dux_new.*ry_new;

% Break up by separation bins
for i = 1:length(dist_axis)
    id = find(dist<dist_bin(i+1) & dist>=dist_bin(i));
    struct_pairs_fast(i).dull = dull_new(id);
    struct_pairs_fast(i).dutt = dutt_new(id);
end
toc
%%
function rx = dist_rx(XI, XJ)
    rx = (XI(:,1) - XJ(:,1)).*cosd(0.5*(XI(:,2) + XJ(:,2)));
end
function ry = dist_ry(XI, XJ)
    ry = (XI(:,2) - XJ(:,2));
end
function du = dist_du(UI, UJ)
    du = UI - UJ;
end
function dist = dist_geo(XI,XJ)
        X = abs(XI(:,1) - XJ(:,1)) .*cosd(0.5*(XI(:,2)+XJ(:,2))) *111321;
        Y = abs(XI(:,2) - XJ(:,2)) *111321;
        dist = sqrt(X.^2 + Y.^2);
end


