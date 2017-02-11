%% Spike and slab no interaction

clear

% cd('c:\Users\Marina Antillon\Google Drive\Dissertation\meta-analysis\spike_slab\spike_slab_LKO')

rng('shuffle')

z=1;

load('groups')

% cd('/Users/Marina/Google Drive/Dissertation/meta-analysis/spike_slab')
data = readtable('../00-data/data_long.txt', 'Delimiter', '\t', 'TreatAsEmpty', 'NA');
% of the validation dataset.

run ../subroutines/prep_data_long.m

% Normalizing covariate values
run ../subroutines/normalizing_data.m

locs_loo = unique(data.location);
locs_loo(strcmp(locs_loo, 'Rehri Goth')) = [];
out = locs_loo(groups_old(:,z),1); % SAVE which location you are taking out each time.
if any(strcmp(out, 'Karachi Slums'))
    out = [out; 'Rehri Goth']; 
end
data(ismember(data.location, out), :) = [];

% Designate a random effects group.
data.region = categorical(data.location);  
data.region = removecats(data.region, cellstr(nominal(out)));
% input.country = double(categorical(datajags.region))';
RE_group_lbl = data.region;
RE_group_num = double(categorical(data.region))';

cd('../03-LKOvalidation')

% Run one jags model.
run ../subroutines/prep_jags_input.m

% Fill in the initial values when the random walk starts with a semi-saturated model.
run ../subroutines/prep_jags_initial_sat.m

% Fill in the initial values when the random walk starts with a null model.
% run ../subroutines/prep_jags_initial_null.m

% Consolidate the initial value structures in one structure

initial(1,1) = initial_sat;
% initial(1,2) = initial_null;
% clearvars -except initial input normalize covars

% Run JAGS 

nburnin = 1e5;
nsamples = 5e5;
nchains = size(initial, 2);
k=z+10; % thread number. In case I am running more than one jags process at a time.
if nchains>1
doparallel = 1; % do not parallelize
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

filename = ['./loo_validation', sprintf('%02.0f', z), '.mat'];

save(filename, 'structArray', 'normalize', 'covar', 'out', 'RE_group_lbl', 'RE_group_num', '-v7.3')

