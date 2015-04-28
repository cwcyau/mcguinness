
% data file
datfile = 'data/Nork.txt'; 

% results output file
resultsfile = 'output/Nork.mat';

% plot of results
plotfile = 'output/Nork.pdf';

% summary file of parameter estimates
summaryfile = 'output/Nork-summary.txt';

% program options
options.maxIters = 5000;
options.burnin = 1000;
options.thin = 1;

options.u = 0.05; % prior mean of large event probability (REMOVE TO USE DEFAULT)
options.v = 0.1^2; % prior variance of large event probability (REMOVE TO USE DEFAULT)

% run bayesian sampling 
options = apsampler(datfile, resultsfile, options);

% make plots
makeplots(datfile, resultsfile, plotfile, summaryfile);





