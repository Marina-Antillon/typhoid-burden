% Post estimation procedure: simulate incidence from the null and the 
% Simulations according to the OLD DATASET 

clear

addpath ../subroutines

rng('shuffle')

% cd('/Users/Marina/Google Drive/Dissertation/meta-analysis/spike_slab')
% cd('c:\Users\Marina Antillon\Google Drive\Dissertation\meta-analysis\spike_slab')

data = readtable('../00-data/data_long.txt', 'Delimiter', '\t', 'TreatAsEmpty', 'NA');

run ../subroutines/prep_data_long.m

% take out adult incidence that was only observed in 15-25 year olds, etc.
% data.cases(data.sel=='CHL Santiago 1984' & data.Age_Group=='Adults 15+') = nan;
% data.cases(data.sel=='CHN Quan 1996' & data.Age_Group=='Adults 15+') = nan;
% data.cases(data.sel=='VNM Hue 2002' & data.Age_Group=='Adults 15+') = nan;

% Designate a random effects group.
% data.region = double(categorical(data.location)); % data.sel; 

% Designate a random effects group.
data.region = data.location; % data.sel; 

% Normalizing covariate values

run ../subroutines/normalizing_data.m

run ../subroutines/prep_data_wide.m

cd('../04-results/')

data_wide.surv = (data_wide.survtype1=='P') ;
% data_wide.surv = zeros(length(data_wide.surveillance1),1);

data_wide.region = double(categorical(data_wide.region)); 

load('../output/jan06_4ages_location_hp_slope.mat') % the results of the model

structArray(1,1).processes=2; % (2) takes into account the obs process

nsamples=2e4; % size(structArray.mu,1);
tic
[prelambda, lambda] = fn_predictions_4ages_slope(data_wide, structArray(1,1), nsamples, covar, 1, 1, 1); 
% Put a 0 in the third to last element to get a prediction for every sample.
toc

save('../output/predictions_model_2proc_jan06RE.mat', 'prelambda', 'lambda')

data_wide.region = zeros(size(data_wide, 1),1);

%% OBSERVED CASES

load('../output/jan06_4ages_location_hp_slope.mat')

% Adjusting for the observation process 
nsamples=1e4;

structArray(1,1).processes=2; % (1) don't take into account the obs process
                              % (2) take into account the obs process.

nsamples=1e4; % size(structArray.mu,1);
tic
[prelambda, lambda] = fn_predictions_4ages_obs(data_wide, structArray(1,1), nsamples, covar, 1, 1, 1); 
% Put a 0 in the third to last element to get a prediction for every sample.
toc

save('../output/predictions_obs_2proc_jan06.mat', 'prelambda', 'lambda')
    
% Without the observation process (from the adjusted model)

structArray(1,1).processes=1; % (1) don't take into account the obs process
                              % (2) take into account the obs process.

nsamples=1e4; % size(structArray.mu,1);
tic
[prelambda, lambda] = fn_predictions_4ages_obs(data_wide, structArray(1,1), nsamples, covar, 1, 1, 1); 
% Put a 0 in the third to last element to get a prediction for every sample.
toc

save('../output/predictions_obs_dproc_only_jan06.mat', 'prelambda', 'lambda')

%% Simulate the distribution of incidence estimated from the NULL MODEL

load('../output/jan06_4ages_location_hp_slope_NULL.mat')

structArray(1,1).processes=2; % (1) don't take into account the obs process
                              % (2) take into account the obs process.

nsamples=1e4; % size(structArray.mu,1);
tic
[prelambda, lambda] = fn_predictions_4ages_slope_NULL(data_wide, structArray(1,1), nsamples, covar, 1, 1, 1); 
% Put a 0 in the third to last element to get a prediction for every sample.
toc

% Save those
save('../output/predictions_NULL_2proc_jan06.mat', 'prelambda', 'lambda')

% I commented it out because I don't need it. But here is the code just in
% case.

% structArray(1,1).processes=1; % (1) don't take into account the obs process
                                % (2) take into account the obs process.
% 
% nsamples=1e4; % size(structArray.mu,1);
% tic
% [prelambda, lambda] = fn_predictions_4ages_slope_NULL(data_wide, structArray(1,1), nsamples, covar, 1, 1, 1); 
% % Put a 0 in the third to last element to get a prediction for every sample.
% toc
% 
% % Save those
% save('../output/predictions_NULL_dproc_only_jan06.mat', 'prelambda', 'lambda')

%% Simulate the distribution of incidence estimated from the ADJUSTED MODEL

load('./output/jan06_4ages_location_hp_slope.mat')

whichrun = 2;
structArray(1,whichrun).processes=2; % (2) takes into account the obs process

nsamples=1e4; % size(structArray.mu,1);
tic
[prelambda, lambda] = fn_predictions_4ages_slope(data_wide, structArray(1,whichrun), nsamples, covar, 1, 1, 1); 
% Put a 0 in the third to last element to get a prediction for every sample.
toc

% Save those
save('../output/predictions_model_2proc_jan06.mat', 'prelambda', 'lambda')

whichrun = 2;
structArray(1,whichrun).processes=1; % (2) takes into account the obs process

nsamples=1e4; % size(structArray.mu,1);
tic
[prelambda, lambda] = fn_predictions_4ages_slope(data_wide, structArray(1,whichrun), nsamples, covar, 1, 1, 1); 
% Put a 0 in the third to last element to get a prediction for every sample.
toc

% Save those
save('../output/predictions_model_dproc_only_jan06.mat', 'prelambda', 'lambda')

