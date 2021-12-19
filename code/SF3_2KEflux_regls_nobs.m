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

%% how we do the inversion.
% Options: LS, NNLS, Reg LS, Reg NNLS
inv_style = 'LS'; 
%lambda = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1] ;
lambda = logspace(-4, -1, 11);
%

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
n=1; 
for m=1:length(lambda)
    
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
       I(j,j) = dk(j)^0.5; 
    end
    
    V = V'; % make a column
    
    Areg = [A; lambda(m)*I];
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


term1(m) = sum((V-Vt./R').^2);
term2(m) = sum((ebs(1:end-1)).^2 .* dk');

%
figure
subplot(211)
semilogx(R, V)
hold all 
semilogx(R, Vt./R')
grid on
title(num2str(lambda(m)))

subplot(212)
semilogx(kf, ebs(1:end-1).*kf')
grid
%savefig(strcat(num2str(lambda(m)), '.pdf'))
print(['lambda_', num2str(lambda(m)), '.png'],'-dpng', '-r400')

end

SpecFlux = flipud(SpecFlux); 


%%
figure
loglog(term1, term2)
text(term1 (1:2:end), term2(1:2:end), num2str(lambda(1:2:end)'))
xlabel('Error norm')
ylabel('Total variation')
print(['L_curve.png'],'-dpng', '-r400')



%%
inv_style = strcat(inv_style, strcat('reg', num2str(lambda))); 

%%
%plots_SF3_fits