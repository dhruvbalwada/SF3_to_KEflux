%% 3rd Order Structure Function to Kinetic Energy Spectral Flux
% Dhruv Balwada, 29 January 2021
% 
% Dhruv Balwada, 6 December 2021
% Improved grid discretization (dk) etc.
%
% Simple least squares

clear all 
close all 

%%

% how we do the inversion.
% Options: LS, NNLS, Reg LS, Reg NNLS
inv_style = 'LS'; 
lambda = 1;
% Experiment name
% Options: LASER, GLAD
experiment = 'LASER'; 

if strcmp(experiment, 'LASER')
    load ../data/LASER_S3_deep500_boot_strap.mat
else
    load ../data/GLAD_S3_deep500_boot_strap.mat
end

s3lll = mean(s3lll,2);
s3ltt = mean(s3ltt,2);

%%

%nsamps = size(s3lll,2);
nsamps=1;

r=dist_axis;
Nr=length(r);

% select part of the r axis that we think has reasonable
% data.

ns=find(dist_axis>=20 ,1);
ne=find(dist_axis<=500e3 ,1,'last');

R=r(ns:ne); 

NR=length(R);
lf=R;

%%
% Define the k axis where we do estimation.

% Same points as r points
%kf=1./lf;

% or 

% Similar points as r points
kf = logspace( log10(1/max(R)), log10(1/min(R)), length(R)-1) ;
dk = diff(kf);
kf = 0.5*(kf(1:end-1) + kf(2:end));

% Code is written in a way that it is assumed that
% kf is decreasing. This matters when doing cumsum
kf = fliplr(kf); 
dk = fliplr(dk);

Nk=length(kf);
%%

% Define matrices
ebs = zeros(Nk+1,nsamps);
S = zeros(Nr, nsamps);
Vt = zeros(NR,nsamps);
SpecFlux = zeros(Nk, nsamps); 

norm_flag=1;

for n=1:nsamps
    
    % SF3 from the data
    S(:,n) =s3lll(:,n)' +s3ltt(:,n)';
    
    % SF3 over the selected range od points
    V=S(ns:ne, n)';
    
    % We will solve Ax=b problem
    % where x is the parameters (energy fluxes), b is SF3
    
    % build the matrix A
    A=zeros(NR,Nk+1);
    for j=1:Nk
        A(:,j)=-4/kf(j)*besselj(1,kf(j)*R)'*dk(j);
    end
    A(:,end)=2*R';
    
    % %% divide Ax=b by r on both sides to remove log dependence
    % normalize the magnitude
    if norm_flag == 1
        for j=1:NR
            A(j,:)=A(j,:)./abs(R(j));
        end
        V=V./abs(R);
    end
    
    % Build the regularization matrices
    I = zeros(Nk+1, Nk+1); 
    O = zeros(Nk+1, 1);
    
    for j=1:Nk
       I(j,j) = dk(j); 
    end
    
    V = V'; % make a column
    
    Areg = [A; lambda*I];
    Vreg = [V; O];
    
    % estimate the epsilons using least squares
    if strcmp(inv_style, 'NNLS')
        ebs(:,n) = lsqnonneg(Areg,Vreg);
    elseif strcmp(inv_style, 'LS')
        ebs(:,n) = Areg\(Vreg);
    end
    
    % What does the reconstructed V look like
    Vt(:,n) =2*ebs(end,n)*R;
    for j=1:Nk
        Vt(:,n)=Vt(:,n)-4*ebs(j,n)./kf(j).*besselj(1,kf(j).*R)'*dk(j);
    end
    
    % Convert V to spectral fluxes   
    SpecFlux(1,n) = -ebs(end,n); 
    
    for j = 2:Nk
        SpecFlux(j,n) = SpecFlux(j-1,n) + ebs(end-j+1, n)*dk(end-j+2);
    end
    
end

SpecFlux = flipud(SpecFlux); 
%% Estimate confidence intervals 

clear CI_ebs CI_Vt

for i = 1:size(ebs,1)
    CI_ebs(:,i) = prctile(ebs(i,:), [99, 1]); 
end

for i = 1:size(Vt,1)
    CI_Vt(:,i) = prctile(Vt(i,:), [99,1]);
end

for i = 1:size(SpecFlux,1)    
    CI_SpecFlux(:,i) = prctile(SpecFlux(i,:), [99,1]);
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
inv_style = strcat(inv_style, strcat('reg', num2str(lambda))); 

%%
plots_SF3_fits