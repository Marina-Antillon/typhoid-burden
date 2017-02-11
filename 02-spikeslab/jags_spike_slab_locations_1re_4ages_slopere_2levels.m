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

% Run one jags model.
data = sortrows(data,{'continent','country'});

run ../subroutines/prep_jags_input.m

input.country = double(categorical(datajags.region))';

input.R = length(unique(categorical(datajags.continent)));
temp = unique([double(categorical(datajags.region)), ...
    double(categorical(datajags.continent))], 'rows');

input.continent = temp(:,2)';

for x=1:length(input.continent)
    input.countrycontinent(x) = sum(input.continent(1:(x-1))==input.continent(x))+1;
end

input.C1 = sum(input.continent==1);
input.C2 = sum(input.continent==2);
input.C3 = sum(input.continent==3);

mustadd = [0, input.C1, input.C1+input.C2];

input.countryposition = input.countrycontinent + mustadd(input.continent);

% Fill in the initial values when the random walk starts with a semi-saturated model.
run ../subroutines/prep_jags_initial_sat_2levels.m
% Fill in the initial values when the random walk starts with a null model.
run ../subroutines/prep_jags_initial_null_2levels.m

% Consolidate the initial value structures in one structure
initial(1,1) = initial_sat;
initial(1,2) = initial_null;
clearvars -except initial input normalize covar

% Run JAGS 
cd('../02-spikeslab/')
nburnin = 1e5;
nsamples = 2.5e5;
nchains = size(initial, 2);
k=round(rand(1)*100); % thread number. In case I am running more than one jags process at a time.
if nchains>1
doparallel = 1; % do not use parallelization
else
doparallel = 0;
end

monitor = {'alpha', 'beta', 'zeta1', 'zeta2', 'zeta3', 'betaprob', ...
    'betau5', 'betao5', 'psi', 'sigma11', 'sigma12', 'sigma13', 'sigma14', ...
    'sigma21', 'sigma22', 'sigma23','sigma24', ...
    'sigma31', 'sigma32', 'sigma33','sigma34', ...
    'omega', 'sigmabeta1', 'sigmabeta2', 'sigmabeta3', 'sigmabeta4', ...
    'betaout1', 'betaout2', 'betaout3', ...
    'gamma', 'gamma_int1',  'gamma_int2', 'gamma_int3', ...
    'gammapick', 'taugamma'}; 

fprintf( 'Running JAGS...\n' );
clock
tic
clear samples stats structArray
[~, ~, structArray] = matjags( ...
    input, ...                     % Observed data   
    fullfile(pwd, './jags_code/tss_int_1re_4ages_hp_slope_taugamma4_2levels.jags'), ...    % File that contains model definition
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

filename = ['./output/' months{today(2)} sprintf('%02.0f', today(3)) '_4ages_location_hp_slope_2levels_continent.mat'];
    
save(filename, 'structArray', 'normalize', 'covar', '-v7.3')

hist([structArray(1,1).omega(:,i), structArray(1,1).betaout1(:,i), structArray(1,1).betaout2(:,i), structArray(1,1).betaout3(:,i)], 100)


regions = {'Global', 'Africa', 'America', 'Asia'};
measures = {'Mean', 'Median', 'Low CI', 'High CI'};

twolevels = nan(16, 4);

for i=1:4
    twolevels((1:4)+(i-1)*4, 1) = mean([structArray(1,1).omega(:,i), structArray(1,1).betaout1(:,i), structArray(1,1).betaout2(:,i), structArray(1,1).betaout3(:,i)])';
    twolevels((1:4)+(i-1)*4, 2) = prctile([structArray(1,1).omega(:,i), structArray(1,1).betaout1(:,i), structArray(1,1).betaout2(:,i), structArray(1,1).betaout3(:,i)], 50)';
    twolevels((1:4)+(i-1)*4, 3) = prctile([structArray(1,1).omega(:,i), structArray(1,1).betaout1(:,i), structArray(1,1).betaout2(:,i), structArray(1,1).betaout3(:,i)], 2.5)';
    twolevels((1:4)+(i-1)*4, 4) = prctile([structArray(1,1).omega(:,i), structArray(1,1).betaout1(:,i), structArray(1,1).betaout2(:,i), structArray(1,1).betaout3(:,i)], 97.5)';
end

twolevels2 = cell(4,4);
twolevels2 = cell2table(twolevels2);

twolevels2.Global = cell(4,1);
twolevels2.Africa = cell(4,1);
twolevels2.America = cell(4,1);
twolevels2.Asia = cell(4,1);

mat2cell(nan(4,4))

for i=1:4
    for j=1:4
        twolevels2(j,i) = cellstr([num2str(twolevels(j+(i-1)*4, 2), '%2.1f'), ' (', num2str(twolevels(j+(i-1)*4, 3), '%2.1f'), ', ', num2str(twolevels(j+(i-1)*4, 4), '%2.1f'), ')']);
    end
end


