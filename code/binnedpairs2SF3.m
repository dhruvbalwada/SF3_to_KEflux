%% Binned pairs to 3rd Order Velocity Structure Functions 
% Dhruv Balwada, 29 January 2021
%
% This script takes the pairs that are separated in different separation
% bins and estimates the third order structure function 
% SF3 = <du_l (du_l^2 + du_t^2)>
%
% Here do bootstrapping to generate confidence intervals on SF3. The N pairs
% present in each bin are resampled to get N random samples (with
% repetition). These resampled data are then used to create an estimate of
% SF3. This process is repeated many times to generate many SF3s.

clear all

%% Load the file where the binned pairs are
load /Users/dhruvbalwada/Work/structure_functions_GOM/structure_pairs.mat

%% Generate many estimates of SF3 with random samples

l=length(dist_axis);

s3lll = zeros(l,2000);
s3llt = zeros(l,2000);

% this loop is extremely slow
for i =1:length(struct_pairs)
    
    npairs = length(struct_pairs(i).dull);
    disp(i)
    for j =1:2000
        y = randsample(npairs, npairs, true);
        
        s3lll(i,j) = nanmean(struct_pairs(i).dull(y).^3);
        s3ltt(i,j) = nanmean(struct_pairs(i).dull(y).*struct_pairs(i).dutt.^2);
    end
end

%% Save to file

save S3_boot_strap.mat s3lll s3ltt dist_axis dist_bin
