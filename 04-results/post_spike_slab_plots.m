%% Plot model incidence

clear

rng('shuffle')

% cd('/Users/Marina/Google Drive/Dissertation/meta-analysis/spike_slab')
% cd('c:\Users\Marina Antillon\Google Drive\Dissertation\meta-analysis\spike_slab')

data = readtable('../00-data/data_long.txt', 'Delimiter', '\t', 'TreatAsEmpty', 'NA');

run ../subroutines/prep_data_long.m

% Designate a random effects group.
data.region = data.location; % data.sel; 

run ../subroutines/prep_data_wide.m

cd('../04-results')

data_wide.surv = (data_wide.survtype1=='P') ;
% data_wide.surv = zeros(length(data_wide.surveillance1),1);
% data_wide.region = double(categorical(data_wide.region)); 
% data_wide.region = zeros(size(data_wide, 1),1);

%% Simulation output for Null, Model, and LKO validation
% Figure 5

load('../output/predictions_NULL_2proc_jan06')

% prelambda(prelambda>1e-1)=1e-1;
% prelambda(prelambda<1e-5)=1e-5;
prelambda_all{1} = prelambda;
% run ../subroutines/plot_IR_violin

load('../output/predictions_model_2proc_jan06')

prelambda = prelambda(:,:,1:1e4);
% prelambda(prelambda>1e-1)=1e-1;
% prelambda(prelambda<1e-5)=1e-5;
prelambda_all{2} = prelambda;
% run ../subroutines/plot_IR_violin

load('../03-LKOvalidation/LOO_sims_slopere_jan06')
prelambda(prelambda>1e-1)=1e-1;
prelambda(prelambda<1e-5)=1e-5;
prelambda_all{3} = prelambda;

prelambda=prelambda_all;

% All results
run ../subroutines/plot_IR_violin_null_model_lko

% set(gca, 'Fontsize', 12, 'FontName','Arial', 'fontWeight', 'Normal')
% saveas(gca,'../figures_dec/three_violins.pdf')

%% Plot with RE
% Figure S3

load('../output/predictions_model_2proc_dec20RE.mat')

% prelambda(prelambda>1e-1)=1e-1;
% prelambda(prelambda<1e-5)=1e-5;
% prelambda_all{1} = prelambda;
run ../subroutines/plot_IR_violin

%% Plot without RE
% Produces Figure 4

clear

addpath ../subroutines

rng('shuffle')

% 
cd('/Users/Marina/Google Drive/Dissertation/meta-analysis/spike_slab')
% cd('c:\Users\Marina Antillon\Google Drive\Dissertation\meta-analysis\spike_slab')

data = readtable('../00-data/data_long.txt', 'Delimiter', '\t', 'TreatAsEmpty', 'NA');

run ../subroutines/prep_data_long.m

% Designate a random effects group.
data.region = data.location; % data.sel; 

% Normalizing covariate values

run ../subroutines/normalizing_data.m

run ../subroutines/prep_data_wide.m

cd('../04-results')

data_wide.surv = (data_wide.survtype1=='P') ;

data_wide.region = zeros(size(data_wide, 1),1);

load('../output/predictions_obs_dproc_only_jan06.mat')
prelambda_obs = prelambda;
lambda_obs = lambda;
load('../output/predictions_model_dproc_only_jan06.mat')

% age groups (second dimension of lambda and prelambda)
% (1) 0-<2 (2) 2-<5 (3) 5-<15 (4) 15+ (5) 0-<5
% (6) 2-<15 (7) 0-<15 (8) 5+ (9) All ages

plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);

age_group_lbl = {'0-<2'; '2-<5'; '5-<15'; '15+'; '0-<5'; ...
                '2-<15'; '0-<15'; '5+'; 'All ages'};
age_order = [1, 5, 2, 7, 6, 3, 8, 4, 9];   
ageval = [1, 2.5, 3.5, 7.5, 8.5, 10, 11, 15, 16]; 

age_group_lbl=age_group_lbl(age_order);

prelambda(prelambda>1e-1) = 1e-1;
prelambda(prelambda<1e-5) = 1e-5;

prelambda_obs(prelambda_obs>1e-1) = 1e-1;
prelambda_obs(prelambda_obs<1e-5) = 1e-5;

prelambda = prelambda(:, age_order, :);
prelambda_obs = prelambda_obs(:, age_order, :);

morethan1 = find(sum(~isnan(prelambda_obs(:, 1:8, 1)), 2)>0);

[~, rankinc] = sort(median(prelambda(morethan1, 9, :), 3), 'descend');
r = 1:length(morethan1);
r(rankinc) = r;

data_wide.location = cellstr(data_wide.location);

data_wide.location(strcmp(data_wide.location, 'Kalkaji')) = {'Delhi'};
data_wide.location(strcmp(data_wide.location, 'Kalamapur')) = {'Dhaka'};
data_wide.location(strcmp(data_wide.location, 'Hechi city')) = {'Hechi City'};
data_wide.location(strcmp(data_wide.location, 'Karachi slums')) = {'Karachi'};
data_wide.location(strcmp(data_wide.location, 'Communes Dong Thap')) = {'Dong Thap'};

mergecountry = table(unique(data_wide.country));
mergecountry.Properties.VariableNames = {'country'} ;

mergecountry.country_name = {'Bangladesh', 'Chile', 'China', 'Egypt', ...
    'Ghana', 'Haiti', 'Indonesia', 'India' 'Kenya', 'Nepal', 'Pakistan', ...
    'Uzbekistan', 'Vietnam', 'South Africa'}';

data_wide = join(data_wide, mergecountry, 'Keys', 'country');

figure('position', [0, 0, 750, 1700], 'color', 'w') 
% Third element is for width, and fourth is for height.
% [ha, pos] = tight_subplot(8,4,[.01 .1],[.1 .1],[.1 .1]);

for i = 1:length(morethan1)

age_group = find(~isnan(prelambda_obs(morethan1(i), :, 1)));
% axes(ha(i))

subplot(4, 8, r(i))
% caxis([0,1])
hold on
plot_variance(ageval(1:8), prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 2.5, 3), ... % 1:length(age_group)
            prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 97.5, 3), [1 0.85 0.85])
plot_variance(ageval(1:8), prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 12.5, 3), ... % 1:length(age_group)
            prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 87.5, 3), [1 0.75 0.75])       
plot_variance(ageval(1:8), prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 22.5, 3), ... % 1:length(age_group)
            prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 77.5, 3), [1 0.65 0.65]) 
plot_variance(ageval(1:8), prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 32.5, 3), ... % 1:length(age_group)
            prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 67.5, 3), [1 0.55 0.55]) 
plot_variance(ageval(1:8), prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 42.5, 3), ... % 1:length(age_group)
            prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 57.5, 3), [1 0.45 0.45]) 
        
plot(ageval(1:8), prctile(1e5*(prelambda(morethan1(i), 1:8, :)), 50, 3), 'Color', [0.75 0 0.25])        
errorbar(ageval(age_group), prctile(1e5*(prelambda_obs(morethan1(i), age_group, :)), 50, 3), ...
        prctile(1e5*(prelambda_obs(morethan1(i), age_group, :)), 50, 3) - prctile(1e5*(prelambda_obs(morethan1(i), age_group, :)), 2.5, 3), ...
        prctile(1e5*(prelambda_obs(morethan1(i), age_group, :)), 97.5, 3) - prctile(1e5*(prelambda_obs(morethan1(i), age_group, :)), 50, 3), ...
        'o', 'Color', [0 0 0], ...
        'MarkerEdgeColor',[0 0 0], ...
        'MarkerFaceColor',[0 0 0], ...
        'MarkerSize', 2, 'LineWidth', 1.5)
    
        
% scatter(log10(prctile(prelambda_obs(i, j, :), 50)), 1:length(age_group), [0.5 0.5 0.5] ,'o', 'filled')
    % to add make proportional to study size, add 'bubblesize' to fourth term
% axis([0 16 -5 0])
set(gca, 'XTick', 0:5:15, 'XTickLabels', {0, 5, 10, 15}', 'XTickLabelRotation', 0, 'fontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 10)
% set(gca, 'YTick', [-4, -3, -2, -1]', 'YTickLabels', {'10', '100', '1,000', '10,000'}', 'YTickLabelRotation', 0, 'fontWeight', 'bold', 'FontName', 'Arial','FontSize', 10)
% set(gca, 'yscale', 'log')

if r(i)<23
    axis([0 16 1 1e4])
else
    axis([0 16 1 1e3])
end

hold off
title({[char(data_wide.country_name(morethan1(i))) , ', ', sprintf('%4.0f', data_wide.midptyr(morethan1(i)))],  ['(', char(data_wide.location(morethan1(i))), ')']}, 'FontName', 'Arial', 'FontSize', 10) % , 'interpreter', 'latex'

end

[ax,h1]=suplabel('Age'); 
[ax,h2]=suplabel('Incidence per 100,000 Person-Years','y', [0.115 0.075 0.8800 0.8800]); 

set(h1,'FontSize',14, 'fontWeight', 'bold','FontName', 'Arial') 
set(h2,'FontSize',14, 'fontWeight', 'bold','FontName', 'Arial') 

% Put a note in figure label that it is adjusting for observation process.

