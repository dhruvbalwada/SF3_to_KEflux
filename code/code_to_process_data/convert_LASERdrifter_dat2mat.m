% Dhruv Balwada
% 3 May 2021 
% Function to extract individual trajectories from the drifter dat files.

clear all
close all

f = fopen('../LASER_SPOT_15min_filtered/laser_spot_drifters_clean_v15.dat');

%%
beep
beep
beep
%%
% There are 4 types of drifter trajs in this data set
% L_ drogued, M_ drogued but cut due to large gap
% U_ undrogued, V_ undrogued but cut due to large gap

nL=0; % drifter counter for different types
nM=0;
nU=0;
nV=0;

while true
    A = fgetl(f);
    
    if A(1) == '%'; continue; end % jump over the header
    
    if ~ischar(A); break; end  %end of file
    
    B = textscan(A, '%s%q%q%f%f%f%f%f%f'); % read line from file as a cell with specified formatting
    
    date_time_str = strcat([B{2}{1},' ', B{3}{1}(1:8)]); % make the time string (throw out milliseconds)
    
    drifter_type = B{1}{1}(1); % read drifter type
    
    switch drifter_type
        case 'L' 
            if nL==0 || ~strcmp(B{1}{1}, drifterL(nL).name)
                nL=nL+1;
                drifterL(nL).name = B{1}{1} ;
                i=1; % counter for the traj point of particular drifter
                disp(drifterL(nL).name)
            end
            
            drifterL(nL).time(i) = datenum(date_time_str);
            drifterL(nL).lat(i)  = B{4};
            drifterL(nL).lon(i)  = B{5};
            drifterL(nL).u(i)    = B{7};
            drifterL(nL).v(i)    = B{8};
            i=i+1;
            
        case 'M'
            if nM==0 || ~strcmp(B{1}{1}, drifterM(nM).name)
                nM=nM+1;
                drifterM(nM).name = B{1}{1} ;
                i=1;
                disp(drifterM(nM).name)
            end
            
            drifterM(nM).time(i) = datenum(date_time_str);
            drifterM(nM).lat(i)  = B{4};
            drifterM(nM).lon(i)  = B{5};
            drifterM(nM).u(i)    = B{7};
            drifterM(nM).v(i)    = B{8};
            i=i+1;
            
        case 'U'
            if nU==0 || ~strcmp(B{1}{1}, drifterU(nU).name)
                nU=nU+1;
                drifterU(nU).name = B{1}{1} ;
                i=1;
                disp(drifterU(nU).name)
            end
            
            drifterU(nU).time(i) = datenum(date_time_str);
            drifterU(nU).lat(i)  = B{4};
            drifterU(nU).lon(i)  = B{5};
            drifterU(nU).u(i)    = B{7};
            drifterU(nU).v(i)    = B{8};
            i=i+1;
                        
        case 'V'
            if nV==0 || ~strcmp(B{1}{1}, drifterV(nV).name)
                nV=nV+1;
                drifterV(nV).name = B{1}{1} ;
                i=1;
                disp(drifterV(nV).name)
            end
            
            drifterV(nV).time(i) = datenum(date_time_str);
            drifterV(nV).lat(i)  = B{4};
            drifterV(nV).lon(i)  = B{5};
            drifterV(nV).u(i)    = B{7};
            drifterV(nV).v(i)    = B{8};
            i=i+1;
                                    
    end
       
end

%% Save different structures in a single mat file

save ../LASER_SPOT_15min_filtered/traj_structs_LASER_15min_03_May_2021.mat drifterL drifterM drifterU drifterV

%%
beep
beep
beep

