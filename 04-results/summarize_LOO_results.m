%% Summarizes LKO validation exercise

% cd('/Users/Marina/Google Drive/Dissertation/meta-analysis/spike_slab')
% cd('c:\Users\Marina Antillon\Google Drive\Dissertation\meta-analysis\spike_slab')

addpath ../subroutines

clear
 
% pull in the randomization table
load('../LKOvalidation/groups')

data = readtable('../00-data/data_long.txt', 'Delimiter', '\t', 'TreatAsEmpty', 'NA');

locs_loo = unique(data.location);
% locs_loo(strcmp(locs_loo, 'Rehri Goth')) = [];

run ../subroutines/prep_data_long.m

cd('../04-results')

% Pull in the MCMC results from each time, rename all stats...

for i=1:7

fname = ['../03-LKOvalidation/loo_validation',  sprintf('%02.0f', i), '.mat' ]; %/spike_slab_LKO
eval(['load ',  fname])
all_posteriors(1, i) = structArray;
norm_stat(1, i) = normalize;
out_all{i,1} = out;

end

logit = @(x) (log(x./(100-x)));

data.wtr_ctr = (logit(data.pipedwater) - normalize.wtr_mean)/normalize.wtr_std;
data.san_ctr = (logit(data.flushtoilets) - normalize.san_mean)/normalize.san_std;
data.gini_ctr = (data.gini_data - normalize.gini_mean)/normalize.gini_std;
data.ext_pov_ctr = (logit(data.pov_hcr_200ppp_perc) - normalize.ext_pov_mean)/normalize.ext_pov_std;
data.road_paved_ctr = (logit(data.road_paved_perc) - normalize.road_paved_mean)/normalize.road_paved_std;
data.eduyrs_f_log_ctr = (log(data.wed) - normalize.eduyrs_f_log_mean)/normalize.eduyrs_f_log_std;
data.stunting_ctr = (logit(data.stunting) - normalize.stunting_mean)/normalize.stunting_std;
data.GDPcap_log_ctr = (log(data.GDPcap) - normalize.GDPcap_log_mean)/normalize.GDPcap_log_std;
data.popdens_log_ctr = (log(data.popdens) - normalize.popdens_log_mean)/normalize.popdens_log_std;
data.floodrisk_log_ctr = (log(data.floodrisk) - normalize.floodrisk_log_mean)/normalize.floodrisk_log_std;
data.water_stress_log_ctr = (log(data.water_blue) - normalize.water_stress_log_mean)/normalize.water_stress_log_std;
data.hivprev_ctr = (log(data.hivprev) - normalize.hivprev_mean)/normalize.hivprev_std;

% Run prep for simulation (make wide data file).

run ../subroutines/prep_data_wide.m

data_wide.surv = (data_wide.survtype1=='P') ;
data_wide.region = zeros(size(data_wide, 1), 1);

% data_wide(data_wide.location=='Transvaal', :) = [];

% Simulate three-four places that were not part of the calibration
% (So 7 rounds of simulations).  

nsamples=1e4; % size(structArray.mu,1);
tic
for i=1:7
out = out_all{i,1}; % locs_loo(groups_old(:,i),1); % SAVE which location you are taking out each time.
    if ismember('Karachi slums', out)
        out = [out; 'Rehri Goth'];
    end 
    
all_posteriors(1, i).processes=2; % (2) takes into account the obs process

[prelambda_tmp{i,1}, lambda_tmp{i,1}] = fn_predictions_4ages_slope(data_wide(ismember(data_wide.location, out), :), all_posteriors(1,i), nsamples, covar, 1, 1, 1); 
% Put a 0 in the third to last element to get a prediction for every sample.
which_data_wide{i,1} = find(ismember(data_wide.location, out));
end
toc
% lambda_tmp stores whole set of nsamples realizations for each place 

% Rearrange it so that the order of the places/studies matches with the
% order in data_wide
lambda = nan(size(data_wide, 1), 9, nsamples); 
for i=1:7
    for j=1:size(which_data_wide{i,1},1)
    lambda(which_data_wide{i,1}(j,1), :, :) = lambda_tmp{i,1}(j, :, :); 
    end
end

prelambda = nan(size(data_wide, 1), 9, nsamples); 
for i=1:7
    for j=1:size(which_data_wide{i,1},1)
    prelambda(which_data_wide{i,1}(j,1), :, :) = prelambda_tmp{i,1}(j, :, :); 
    end
end

save('LOO_sims_slopere_jan06.mat', 'prelambda', 'lambda') 

%% Table of the covariates chosen by spike-and-slab
% Makes Figures S5 and S6

for j=1:7
    for k = 1:3
        forgraph(length(covar):-1:1, j, k) = mean(all_posteriors(1,j).gammapick(1:5e5, :)==k)';
    end
end

COVAR = {'Population Density', 'GDP per Capita', 'Gini Coefficient', ...
        'Access to Piped Water', 'Access to Flush Toilets', 'Years of Education, Women', 'Percent Paved Roads', ...
        'People Living in Extreme Poverty', 'Prevalence of Stunting', 'Flood Risk', 'Prevalence of HIV',...
        'Water Stress'}; %
    
% figure('position', [0, 0, 1000, 450], 'color', 'w')
% [ha, pos] = tight_subplot(1,2,[.01 .01],[.1 .1],[.075 .075]); 
% 
% axes(ha(1))

% subplot(1,2,1)
plotBarStackGroups(forgraph, COVAR(end:-1:1), 'horizontal')
set(gca, 'Fontsize', 12)
xlabel('Proportion of models', 'FontSize', 12)
legend({'Exclude', 'Intercept', 'Intercept & Slopes'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12) 
colormap([0.3, 0.3, 0.3; 0.6, 0.6, 0.6])

set(gca, 'Fontsize', 12, 'FontName','Arial', 'fontWeight', 'Bold')
% saveas(gca,'../figure_dec/stacked_selection_LKO.pdf')
% heatmap: how often are two variables present together?
% Is a better question: of all the times that a variable P is present, what
% percent of the time is the other variable Q also present?

% How many variables does each model have? Most models are not the null model.
colormap(gray(7))
% subplot(1,2,2)
h = hist([sum(all_posteriors(1,1).gammapick(1:5e5, :)>1, 2), ...
    sum(all_posteriors(1,2).gammapick(1:5e5, :)>1, 2) ...
    sum(all_posteriors(1,3).gammapick(1:5e5, :)>1, 2) ...
    sum(all_posteriors(1,4).gammapick(1:5e5, :)>1, 2) ...
    sum(all_posteriors(1,5).gammapick(1:5e5, :)>1, 2) ...
    sum(all_posteriors(1,6).gammapick(1:5e5, :)>1, 2) ...
    sum(all_posteriors(1,7).gammapick(1:5e5, :)>1, 2)], 0:12); 
% title('Number of predictors included in each model', 'Fontsize', 14)

bar(0:12, h/length(structArray(1,1).gammapick))
xlabel('Number of Predictors', 'Fontsize', 12)
ylabel('Proportion of Models', 'Fontsize', 12) 

legend({'Leave-k-out group #1', 'Leave-k-out group #2', 'Leave-k-out group #3', ...
    'Leave-k-out group #4', 'Leave-k-out group #5', 'Leave-k-out group #6', 'Leave-k-out group #7'}, ...
    'Location', 'NorthEast', 'Orientation', 'Vertical', 'Fontsize', 12)

set(gca, 'Fontsize', 12, 'FontName','Arial', 'fontWeight', 'Bold')


