## General structure of code

The general workflow here is: 

**Access data**:
- Download the data from the Carthe website: https://data.gulfresearchinitiative.org/ Data can be found by searching for GLAD and LASER drifters.  
- Conver the dat files to mat files. Code to do this is present in `code_to_process_data`. 

**Key calculations**:
- Calculate the ocean depth at each points of the trajectory and add to data file. Done in `trajectory_depth.m`. This is used in the next step to only keep data that is in deep ocean. 
- Convert trajectories to binned pairs. (Use the vectorized version) `trajectories2binnedpairs_vectorized.m`. Generates files with names like `structure_pairs*.mat`.
    - In ths past this was done using the script called `trajectories2binnedpairs.m`, which used loops and relied on a function called `calculate_seperation_timeseries.m`. This is very slow. 
- Calculate structure functions (2nd and 3rd order) from the binned pairs, and store the mean values in `SF_*.mat`. Dones in `pairtimes2SF.m`
- Do the helmholtz decomposition of SF2. At the moment done using a function called `helmholtz_decompose` in `plot_SFs.m`. SF2 decomposed not saved anywhere.
- Generate many SF3 samples using bootstrapping. Done in `binnedpairs2bootstrappedSF3.m` (the copy of this script was created to be able to run in parallel since this code takes a long time.)
- Convert from SF3 to KE flux. Done in `SF3_2KEflux_ls.m`.
- 

**Plotting functions**:
- trajectories_movie.m, trajectories_plot.m, LASER_trajectories_movie.m: plot trajectories on a map. 
- plot_SFs.m: plots of SF2, the data distribution, pdfs of du, etc. 
- SF3_2KEfluxm: plots of SF3 and associated calculations.
