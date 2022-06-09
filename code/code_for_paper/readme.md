List of what different functions do:

Important scripts: \
1. trajectories2binnedpairs_vectorized.m : Converts matrices of trajectories into all pairs at each time. Matrices of trajectories has time of experiment along one axis and drifter number along the other axis. So if some drifters starts later, it will have some nans initially. The code to generate these are present in the code_to_process_data folder. \
\
2. pairtime2SF.m : Go from pairs as a function of time to pairs in each separation bins. Then the moments are calculated too, but we want to do that with bootstrapping to get error bars (that is in errorbarsSF*).\
\
3. SF3_2KEflux_ls.m : Main function that converts from SF3 to KE flux. It takes in bootstrapped samples. Also choice of method to do the inversion (least squares, least squares with non-negative, with regularization etc can be picked). \
\
4. helmholtz_decompose.m -  Function to do the Helmholtz decomposition.\
\
5. plots_SF3_fits_*.m : Code to plot panels for SF3.\
\
6. plots_SF2s.m : Code to plots SF2.\
\
7. plots_datadist.m : Plot how the data is distributed.\
\
8. lagrangian_rotary_spectra.m : Compute the frequency spectra from the data.\
\
9. errorbarsSF*.m : Do modified block bootstrapping to estimate errorbars and spread using resampling. \
\
Misc scripts:\
1. trajectory_depth.m : Estimate the depth of the ocean at each lat and long where the drifter is.\
2. trajectories_plot.m: \
3. trajectories_movie.m\
4. trajectories2binnedpairs.m - old unvectorized version, that works very slowly. \
5. test_new_SFcalc.m - just a simple test to see if the vecotrized way of doing the above function works.\
6. SF3_2KEflux_regls_nobs.m, SF3_2KEflux_regls.m  - old versions of code to convert from SF3 to KE flux.\
7. SF2_and_helmholtz_decompos_old.m - Really old code to compute structure functions and make plots. Don\'92t use this.\
8. calculate_separation_timeseries.m - Really old code to compute separation time series. Do not use. Use vectorized method instead.}
