%20190911
% Downloaded from:https://github.com/hannnwang/Anisotropic-Helmholtz-decomposition-for-Lagrangian-data/blob/master/LASER/dataformating_all.m
% on 3 May 2021
%Format the recordings from all launches in LASER drifters.
%To run this code, you need to have laserspotdrifterscleanv15_all.mat ready, the directions 
% to which is written in line 11. You can also write your own version if you find it confusing.
%However you do it, in the end, there should be four matrices trajmat_X, trajmat_Y, trajmat_U, 
% trajmat_V that all have dimensions ntime*ndrifters, where ntime is the number of snapshots 
% available in LASER, and ndrifters is the number of drifters available in LASER. Abscent 
%recordings are denoted by NaN. The matrix trajmat_X denotes the recorded longitude, trajmat_Y 
%denots the recorded latitude, trajmat_U denotes the recorded zonal velocities, and trajmat_V 
% denotes the recorded meriodinal velocities. Importantly, the time column is sorted in ascending 
% order. All units are the same as in https://data.gulfresearchinitiative.org/data/R4.x265.237:0001
% The four matrices are then saved as LASER_all.mat, which will be used in the decomposition codes.


clear all
close all

load ../LASER_SPOT_15min_filtered/laserspotdrifterscleanv15_all.mat 
%This .mat file is converted directly from laserspotdrifterscleanv15.dat downloaded 
% from https://data.gulfresearchinitiative.org/data/R4.x265.237:0001
%To convert laserspotdrifterscleanv15.dat into laserspotdrifterscleanv15_all.mat, 
% simply run Matlab, select "Import Data", and load laserspotdrifterscleanv15.dat. 
% When you see the table in Matlab, select the "Range" as A6:BB3564200 (which should 
% be the default Matlab selection), "Output Type" as "String array" (which is likely not 
% the Matlab default option), and import. Rename the imported matrix as "laserspotdrifterscleanv15",
% and save as laserspotdrifterscleanv15_all.mat.

%%
laserspotdrifterscleanv15(:,1)=erase(laserspotdrifterscleanv15(:,1),"L_");

%%

ID_all=unique(laserspotdrifterscleanv15(:,1));
% %Find the ones from P1 launch only
% ID_ifP1=ismember(ID_all,ID_P1);
% 
% laserP1=laserspotdrifterscleanv15(ID_ifP1,:);

%%
%Simplify the IDs
IDsimple=1:length(ID_all);

%laser_all_ordered=changem(laserspotdrifterscleanv15,IDsimple,(ID_all));
laser_all_ordered = laserspotdrifterscleanv15;
% %Sort according to time
% %sort from date first:

% laserP1=sortrows(laserP1,2);
% %then, sort from time (tricker)
% for 
%%

%Sort according to time
%find month, day, hour and minute
date=char(laser_all_ordered(:,2));
month=str2num(date(:,6:7));
day=str2num(date(:,end-1:end));
clearvars date;
time=char(laser_all_ordered(:,3));
hour=str2num(time(:,1:2));
minute=str2num(time(:,4:5));
clearvars time;
for i=1:length(month)
    if month(i)==1
        timearray(i)=minute(i)+60*hour(i)+60*24*day(i)+...
    60*24*0;
    elseif month(i)==2
        timearray(i)=minute(i)+60*hour(i)+60*24*day(i)+...
    60*24*31;
    elseif month(i)==3
        timearray(i)=minute(i)+60*hour(i)+60*24*day(i)+...
    60*24*(31+29);
  elseif month(i)==4
        timearray(i)=minute(i)+60*hour(i)+60*24*day(i)+...
    60*24*(31+29+31);
    else
        warning('something is wrong in month!')
    end
end
timearray=(timearray-min(timearray))/15+1;

%%
laser_all(:,1)=str2num(char(laser_all_ordered(:,1)));%floater ID
laser_all(:,2)=timearray;%Time array
laser_all(:,3)=str2num(char(laser_all_ordered(:,4)));%Lat
laser_all(:,4)=str2num(char(laser_all_ordered(:,5)));%lon
laser_all(:,5)=str2num(char(laser_all_ordered(:,7)));%U
laser_all(:,6)=str2num(char(laser_all_ordered(:,8)));%V

%laser_all=sortrows(laser_all,2);
%%
nfl=length(ID_all);
%Save as a single huge variable, like Dhruv did.
%(using matrix instead of struct array)
trajmat_X=nan*zeros(max(timearray),nfl);
trajmat_Y=trajmat_X;
trajmat_U=trajmat_X;
trajmat_V=trajmat_X;
for i=1:length(day)
    iflt=laser_all(i,1);
    it=laser_all(i,2);
    trajmat_X(it,iflt)=laser_all(i,4);%Note lat and lon
    trajmat_Y(it,iflt)=laser_all(i,3);
    trajmat_U(it,iflt)=laser_all(i,5);
    trajmat_V(it,iflt)=laser_all(i,6);
end

%%
save('LASER_all_allvariables.mat')
clearvars laserspotdrifterscleanv15 laser_all_ordered day hour ID_all IDsimple ...
    minute month timearray
save('LASER_all.mat')

%%
