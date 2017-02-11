% Post-estimation procedure: Summarize MCMC (JAGS) output.

% Calculates the WAIC (if we have several models) and convergence diagnostics
% Clustered/stacked bar graph of covariates selection for any number of chains
% Distribution of the number of covariates that are included in each model.
% Summarizes the gamma and gamma_int# values (conditional on inclusion)
% Heatmap of how often any two variables are present together in the model.
cd('/Users/Marina/Google Drive/Dissertation/meta-analysis/spike_slab')

load('./output/jan06_4ages_location_hp_slope.mat')

%% Clustered, stacked bar graph of covariate selection

for j=1:2
    for k = 1:3
        forgraph(length(covar):-1:1, j, k) = mean(structArray(1,j).gammapick==k)';
    end
end

COVAR = {'Population Density', 'GDP per Capita', 'Gini Coefficient', ...
        'Access to Piped Water', 'Access to Flush Toilets', 'Years of Education, Women', 'Percent Paved Roads', ...
        'People Living in Extreme Poverty', 'Prevalence of Stunting', 'Flood Risk', 'Prevalence of HIV',...
        'Water Stress'}; %

subplot(3, 1, 1)

plotBarStackGroups(forgraph, flip(COVAR), 'horizontal')
% set(gca, 'XTickLabel', COVAR, 'XTickLabelRotation', 45, 'Fontsize', 12)
colormap([0, 0, 0; 0.3, 0.3, 0.3; 0.75, 0.75, 0.75])
xlabel('Proportion of models', 'Fontsize', 12)
legend({'Exclude', 'Intercept', 'Intercept & Slopes'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Fontsize', 12) 
set(gca, 'Fontsize', 12, 'FontName','Arial', 'fontWeight', 'Bold')
set(gca, 'XTick', 0:0.2:1, 'XTickLabel', {'0', '0.20', '0.40', '0.60', '0.80', '1.00'}, 'FontSize', 12)
text(-0.7, 13, 'A', 'FontSize', 12, 'FontWeight', 'bold', 'FontName','Arial') 

% Plot age IRR's

subplot(3, 1, 2)

colormap([0.3, 0.3, 0.3; 0.6, 0.6, 0.6])
h = hist([sum(structArray(1,1).gammapick>1, 2), sum(structArray(1,2).gammapick>1, 2)],0:1:12);

bar(h/length(structArray(1,1).gammapick))
set(gca,'XTick',1:13)

set(gca,'XTickLabel',{0:12})
xlabel('Number of Predictors')
ylabel('Proportion of Models') 
legend({'Chain 1', 'Chain 2'}, 'Location', 'Northeast', 'Orientation', 'Vertical', 'Fontsize', 12)
set(gca, 'Fontsize', 12, 'FontName', 'Arial', 'fontWeight', 'Bold')
set(gca, 'YTick', 0:0.05:0.3, 'YTickLabel', {'0', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30'}, 'FontSize', 12)

text(-1.5, 0.25, 'B', 'FontSize', 12, 'FontWeight', 'bold', 'FontName','Arial') 

subplot(3, 1, 3)

violincolor = [0.5, 0.5, 0.5];
markercolor = [0.5, 0.5, 0.5];

violin(log10(exp(structArray(1,1).beta)), 'facecolor', violincolor, 'bw', 0.3, ...
    'edgecolor', [], 'facealpha', 0.3, 'scale', 1/2, 'mc', [], 'medc', []);
hline = refline(0, 0)
hline.Color = [0.5 0.5 0.5];
hline.LineStyle = '--';
set(gca, 'XTick', 1:3, 'XTickLabel', {'Ages <2'; 'Ages 2-4'; 'Ages ?15'}, 'TickLength',[0.01 0.01], 'fontWeight', 'bold', 'XTickLabelRotation', 0)
axis([0.5 3.5 -3 1])
set(gca, 'YTick', -3:1, 'YTickLabel', [0.001; 0.01; 0.1; 1; 10], 'FontSize', 12)
% title({'Age-specific', 'Incidence Rate Ratios'} , 'FontSize', 10)
ylabel({'Incidence Rate Ratio'}, 'fontWeight', 'bold', 'FontSize', 12)

text(-0.0001, 1, 'C', 'FontSize', 12, 'FontWeight', 'bold', 'FontName','Arial') 

%% How often is the variable chosen?

prctile(sum(structArray(1,1).gammapick>1, 2), [2.5 50 97.5])
prctile(sum(structArray(1,2).gammapick>1, 2), [2.5 50 97.5])

%% Violin plots to show the predictor coefficients (when they are in the model).
% This is Figure S2.

figure('position', [0, 0, 1000, 1500], 'color', 'w')
[ha, pos] = tight_subplot(4,1,[.01 .01],[.25 .01],[.075 .075]); 

violincolor = [0.5, 0.5, 0.5];
markercolor = [0.75, 0 , 0];

% Intercept
axes(ha(4));

for j=1:12
post_gamma{j,1} = structArray(1,1).gamma(structArray(1,1).gammapick(:,j)>1,j);
post_gamma_sum(:,j) = post_gamma{j,1}(randsample(1:length(post_gamma{j,1}), 1000))';
end
bubblesize = 100*mean(structArray(1,1).gammapick>1)';

violin(post_gamma_sum,'x', 1:12, 'mc', [], 'medc', [], 'bw', 0.2, ...
'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', .4); 

hold on
scatter(1:12, median(post_gamma_sum,1)', bubblesize, 'filled', ...
'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)

refline(0,0)
hold off

% Slope 1
axes(ha(3));

for j=1:12
post_gamma{j,1} = structArray(1,1).gamma_int1(structArray(1,1).gammapick(:,j)==3,j);
post_gamma_sum(:,j) = post_gamma{j,1}(randsample(1:length(post_gamma{j,1}), 1000))';
end
bubblesize = 100*mean(structArray(1,1).gammapick==3)';

violin(post_gamma_sum,'x', 1:12, 'mc', [], 'medc', [], 'bw', 0.2, ...
'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', .4); 

hold on
scatter(1:12, median(post_gamma_sum,1)', bubblesize, 'filled', ...
'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)

refline(0,0)
hold off

% Slope 2
axes(ha(2));

for j=1:12
post_gamma{j,1} = structArray(1,1).gamma_int2(structArray(1,1).gammapick(:,j)==3,j);
post_gamma_sum(:,j) = post_gamma{j,1}(randsample(1:length(post_gamma{j,1}), 1000))';
end
bubblesize = 100*mean(structArray(1,1).gammapick==3)';

violin(post_gamma_sum,'x', 1:12, 'mc', [], 'medc', [], 'bw', 0.2, ...
'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', .4); 

hold on
scatter(1:12, median(post_gamma_sum,1)', bubblesize, 'filled', ...
'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)

refline(0,0)
hold off

% Slope 3
axes(ha(1));

for j=1:12
post_gamma{j,1} = structArray(1,1).gamma_int3(structArray(1,1).gammapick(:,j)==3,j);
post_gamma_sum(:,j) = post_gamma{j,1}(randsample(1:length(post_gamma{j,1}), 1000))';
end
bubblesize = 100*mean(structArray(1,1).gammapick==3)';

violin(post_gamma_sum,'x', 1:12, 'mc', [], 'medc', [], 'bw', 0.2, ...
'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', .4); 

hold on
scatter(1:12, median(post_gamma_sum, 1)', bubblesize, 'filled', ...
'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)

refline(0,0)
hold off

COVAR = {'Population Density', 'GDP per Capita', 'Gini Coefficient', ...
        'Access to Piped Water', 'Access to Flush Toilets', 'Years of Education, Women', 'Percent Paved Roads', ...
        'People Living in Extreme Poverty', 'Prevalence of Stunting', 'Flood Risk', 'Prevalence of HIV',...
        'Water Stress'}; %
    
set(ha(1:4),'XTickLabel',''); set(ha(1:4),'YTickLabel','')
set(ha(1:4),'Xtick', 1:12); set(ha(1:4),'Ytick', [-5 0 5])
set(ha(1:4),'XLim', [0 13]); set(ha(1:4),'YLim', [-6 6])

set(ha(4),'XTickLabel',  COVAR, 'XTickLabelRotation', 90)
set(ha(1:4),'YTickLabel', {'-5', '0', '5'}, 'YTickLabelRotation', 90)

ylabel(ha(4), 'Intercept')
ylabel(ha(3), 'Slope 1')
ylabel(ha(2), 'Slope 2')
ylabel(ha(1), 'Slope 3')

% Only label the variables in the bottom one. Rotate them by 90 degrees
% Label the y-axis and rotate the numbers by 90 degrees

%% Sigma covariance

sigma_mat(1:4, 1) = median(structArray(1,1).sigma1, 1)';
sigma_mat(1:4, 2) = median(structArray(1,1).sigma2, 1)';
sigma_mat(1:4, 3) = median(structArray(1,1).sigma3, 1)';
sigma_mat(1:4, 4) = median(structArray(1,1).sigma4, 1)';
cgo = clustergram(sigma_mat) % , 'ColorMap', 'redbluecmap'
HeatMap(sigma_mat)

mu = prctile(structArray(1,1).mu, [50 2.5 97.5]);
beta1 = prctile(structArray(1,1).beta(:,1), [50 2.5 97.5]);
beta2 = prctile(structArray(1,1).beta(:,2), [50 2.5 97.5]);
beta3 = prctile(structArray(1,1).beta(:,3), [50 2.5 97.5]);

[mean(structArray(1,1).sigma1, 1); mean(structArray(1,1).sigma2, 1); 
    mean(structArray(1,1).sigma3, 1); mean(structArray(1,1).sigma4, 1)]

