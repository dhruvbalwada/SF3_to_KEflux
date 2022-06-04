clear all
close all

f = fopen('../GLAD_15min_filtered/GLAD_15min_filtered.dat');

%%
beep 
beep
beep
%%
n=0; % drifter counter

while true
    A = fgetl(f);
    
    if A(1) == '%'; continue; end % jump over the header
            
    if ~ischar(A); break; end  %end of file
    
     
    B = textscan(A, '%s%q%q%f%f%f%f%f%f') ;
        
    m = strcat([B{2}{1},' ', B{3}{1}(1:8)]);
    
    if n==0 || ~strcmp(B{1}{1}, drifter(n).name)
        n=n+1;
        drifter(n).name = B{1}{1} ;
        i=1;
        disp(n)
    end
    

    drifter(n).time(i) = datenum(m);
    drifter(n).lat(i) = B{4};
    drifter(n).lon(i) = B{5};
    drifter(n).u(i) = B{7};
    drifter(n).v(i) = B{8};
    i=i+1;
    
end

%% 

save traj_struct_GLAD_15min_03_May_2021.mat drifter

%%
beep 
beep
beep

