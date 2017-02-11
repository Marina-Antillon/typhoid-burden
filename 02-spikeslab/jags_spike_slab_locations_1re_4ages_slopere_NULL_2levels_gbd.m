%% Spike and slab no interaction

clear

rng('shuffle')

% cd('/Users/Marina/Google Drive/Dissertation/meta-analysis/spike_slab')
% cd('c:\Users\Marina Antillon\Google Drive\Dissertation\meta-analysis\spike_slab')
data = readtable('../00-data/data_long.txt', 'Delimiter', '\t', 'TreatAsEmpty', 'NA');

run ../subroutines/prep_data_long.m

% Designate a random effects group.
data.region = data.location; % data.sel; 

% Normalizing covariate values

run ../subroutines/normalizing_data.m
cd('../02-spikeslab')

% Run one jags model.
data = sortrows(data,{'GBD_super_regions','country'});

run ../subroutines/prep_jags_input.m

input.country = double(categorical(datajags.region))';

input.R = length(unique(categorical(datajags.GBD_super_regions)));
temp = unique([double(categorical(datajags.region)), ...
    double(categorical(datajags.GBD_super_regions))], 'rows');

input.GBD_super_regions = temp(:,2)';

for x=1:length(input.GBD_super_regions)
    input.countrycontinent(x) = sum(input.GBD_super_regions(1:(x-1))==input.GBD_super_regions(x))+1;
end

input.C1 = sum(input.GBD_super_regions==1);
input.C2 = sum(input.GBD_super_regions==2);
input.C3 = sum(input.GBD_super_regions==3);
input.C4 = sum(input.GBD_super_regions==4);
input.C5 = sum(input.GBD_super_regions==5);
input.C6 = sum(input.GBD_super_regions==6);

mustadd = [0, input.C1, input.C1+input.C2, input.C1+input.C2+input.C3, ...
    input.C1+input.C2+input.C3+input.C4, input.C1+input.C2+input.C3+input.C4+input.C5];

input.countryposition = input.countrycontinent + mustadd(input.GBD_super_regions);

% Fill in the initial values when the random walk starts with a null model.
run ../subroutines/prep_jags_initial_null_2levels2.m

% Consolidate the initial value structures in one structure
initial(1,1) = initial_null;
initial = rmfield(initial, {'gammapick', 'gammaprior1', 'gammaprior2', 'gammaprior3', 'gammaprior4'});
clearvars -except initial input normalize covars

% Run JAGS 
nburnin = 1e4;
nsamples = 1e5;
nchains = size(initial, 2);
k=round(randsample(1:100,1)); % thread number. In case I am running more than one jags process at a time.
if nchains>1
doparallel = 1; % do not use parallelization
else
doparallel = 0;
end
monitor = {'alpha', 'beta', 'zeta1', 'zeta2', 'zeta3', 'betaprob', ...
    'betau5', 'betao5', 'psi', 'sigma11', 'sigma12', 'sigma13', 'sigma14', ...
    'sigma21', 'sigma22', 'sigma23','sigma24', ...
    'sigma31', 'sigma32', 'sigma33','sigma34', ...
    'sigma41', 'sigma42', 'sigma43','sigma44', ...
    'sigma51', 'sigma52', 'sigma53','sigma54', ...
    'sigma61', 'sigma62', 'sigma63','sigma64', ...
    'omega', ...
    'sigmabeta1', 'sigmabeta2', 'sigmabeta3', 'sigmabeta4', ...
    'betaout1', 'betaout2', 'betaout3', 'betaout4', 'betaout5', 'betaout6'}; 

fprintf( 'Running JAGS...\n' );
clock
tic
clear samples stats structArray
[~, ~, structArray] = matjags( ...
    input, ...                     % Observed data   
    fullfile(pwd, './jags_code/tss_int_1re_4ages_hp_slope_NULL_2levels_gdb.jags'), ...    % File that contains model definition
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

filename = ['./output/' months{today(2)} sprintf('%2.0f', today(3)) '_4ages_location_hp_slope_NULL_2levels2.mat'];
    
save(filename, 'structArray', '-v7.3')

%% Make a couple of plots to show that the priors for beta (incidence, and age IRR's) are not significant

hist([structArray.omega(:,i), structArray.betaout1(:,i), structArray.betaout2(:,i), structArray.betaout3(:,i),...
    structArray.betaout4(:,i), structArray.betaout5(:,i), structArray.betaout6(:,i),], 100)

regions = {'Global', 'Africa', 'America', 'Asia'};
measures = {'Mean', 'Median', 'Low CI', 'High CI'};

twolevels = nan(28, 4);

for i=1:4
    twolevels((1:7)+(i-1)*7, 1) = mean([structArray.omega(:,i), structArray.betaout1(:,i), structArray.betaout2(:,i), structArray.betaout3(:,i), ...
         structArray.betaout4(:,i), structArray.betaout5(:,i), structArray.betaout6(:,i)])';
    twolevels((1:7)+(i-1)*7, 2) = prctile([structArray.omega(:,i), structArray.betaout1(:,i), structArray.betaout2(:,i), structArray.betaout3(:,i), ...
                 structArray.betaout4(:,i), structArray.betaout5(:,i), structArray.betaout6(:,i)], 50)';
    twolevels((1:7)+(i-1)*7, 3) = prctile([structArray.omega(:,i), structArray.betaout1(:,i), structArray.betaout2(:,i), structArray.betaout3(:,i), ...
                 structArray.betaout4(:,i), structArray.betaout5(:,i), structArray.betaout6(:,i)], 2.5)';
    twolevels((1:7)+(i-1)*7, 4) = prctile([structArray.omega(:,i), structArray.betaout1(:,i), structArray.betaout2(:,i), structArray.betaout3(:,i), ...
                 structArray.betaout4(:,i), structArray.betaout5(:,i), structArray.betaout6(:,i)], 97.5)';
end

twolevels2 = cell(4,4);
twolevels2 = cell2table(twolevels2);

mat2cell(nan(4,4))

for i=1:4
    for j=1:4
        twolevels2(j,i) = cellstr([num2str(twolevels(j+(i-1)*4, 2), '%2.1f'), ' (', num2str(twolevels(j+(i-1)*4, 3), '%2.1f'), ', ', num2str(twolevels(j+(i-1)*4, 4), '%2.1f'), ')']);
    end
end


