
% data file
datfile = 'data/Nork.txt'; 

% results output file
resultsfile = 'output/Nork.mat';

% plot of results
plotfile = 'output/Nork.pdf';

% program options
options.maxIters = 50000;
options.burnin = 1000;
options.thin = 1;

% run bayesian sampling 
apsampler(datfile, resultsfile, options);

% make plots
makeplots(datfile, resultsfile, plotfile);





