function [sep] = calculate_seperation_timeseries(traj)

nflts = size(traj.X,2);
ndays = size(traj.X,1);

n=1;
% npairs = factorial(nflts)/factorial(nflts-2)/factorial(2);
npairs = nflts*(nflts-1)/2; 

sep(npairs).X = nan*ones(ndays,1);
sep(npairs).Y = nan*ones(ndays,1);
sep(npairs).dist = nan*ones(ndays,1);

for i=1:nflts-1
    for j = i+1:nflts
        sep(n).X = abs(traj.X(:,i) - traj.X(:,j)).*cosd(0.5*(traj.Y(:,i)+traj.Y(:,j)))*111321;
        sep(n).Y = abs(traj.Y(:,i) - traj.Y(:,j))*111321;
        sep(n).dist = sqrt(sep(n).X.^2 + sep(n).Y.^2);
        sep(n).X1 = traj.X(:,i);
        sep(n).X2 = traj.X(:,j);
        sep(n).Y1 = traj.Y(:,i);
        sep(n).Y2 = traj.Y(:,j);
        sep(n).U1 = traj.U(:,i);
        sep(n).U2 = traj.U(:,j);
        sep(n).V1 = traj.V(:,i);
        sep(n).V2 = traj.V(:,j);
        
        % If we want acceleration SF
        %sep(n).Au1 = traj.Au(:,i);
        %sep(n).Au2 = traj.Au(:,j);
        %sep(n).Av1 = traj.Av(:,i);
        %sep(n).Av2 = traj.Av(:,j);
        
        n = n+1;
    end
end
