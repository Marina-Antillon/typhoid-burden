%% Spike and slab no interaction

clear

rng('shuffle')

% cd('/Users/Marina/Desktop/typhoid-burden/02-spikeslab')

data = readtable('../00-data/data_long.txt', 'Delimiter', '\t', 'TreatAsEmpty', 'NA');

run ../subroutines/prep_data_long.m

% Designate a random effects group.

data.region = data.location; % data.sel; 

% Normalizing covariate values

run ../subroutines/normalizing_data.m

% Run one jags model.

run ../subroutines/prep_jags_input.m

% Fill in the initial values when the random walk starts with a semi-saturated model.

run ../subroutines/prep_jags_initial_sat.m

% Fill in the initial values when the random walk starts with a null model.

run ../subroutines/prep_jags_initial_null.m

% Consolidate the initial value structures in one structure

initial(1,1) = initial_sat;
initial(1,2) = initial_null;
clearvars -except initial input normalize covar

% Run JAGS 
cd('../02-spikeslab/')
nburnin = 1e5;
nsamples = 1e6;
nchains = size(initial, 2);
k=round(rand(1)*100); % thread number. In case I am running more than one jags process at a time.
if nchains>1
doparallel = 1; % do not use parallelization
else
doparallel = 0;
end
monitor = {'mu', 'alpha', 'beta', 'zeta1', 'zeta2', 'zeta3', 'betaprob', ...
    'betau5', 'betao5', 'psi', 'sigma1', 'sigma2', 'sigma3', 'sigma4',...
    'gamma', 'gamma_int1',  'gamma_int2', 'gamma_int3', ...
    'gammapick', 'taugamma'}; % 'mugamma',

fprintf( 'Running JAGS...\n' );
clock
tic
clear samples stats structArray
[~, ~, structArray] = matjags( ...
    input, ...                     % Observed data   
    fullfile(pwd, './jags_code/tss_int_1re_4ages_hp_slope_taugamma4.jags'), ...    % File that contains model definition
    initial, k, ...                     % Initial values for latent variables
    'doparallel' , doparallel, ...      % Parallelization flag
    'nchains', nchains,...              % Number of MCMC chains
    'nburnin', nburnin,...              % Number of burnin steps
    'nsamples', nsamples, ...           % Number of samples to extract
    'thin', 10, ...                      % Thinning parameter
    'dic', 1, ...                       % Do the DIC?
    'monitorparams', monitor, ...       % List of latent variables to monitor
    'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
    'verbosity' , 2 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 0 );                    % clean up of temporary files?
toc

% Make a name for the file
months = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'}; 
today = clock;

filename = ['./output/' months{today(2)} sprintf('%2.0f', today(3)) '_4ages_location_hp_slope.mat'];
    
save(filename, 'structArray', 'normalize', 'covar', '-v7.3')

