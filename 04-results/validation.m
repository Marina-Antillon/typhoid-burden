%% (In)Validation of burden model using TSAP data

% clear

rng('shuffle')

addpath ../subroutines

cd('/Users/Marina/Desktop/typhoid-burden/04-results')

load('tsap_logincidence')
load('tsap_simulations')

prelambda(prelambda<1e-5)=1e-5;
prelambda(prelambda>1e-1)=1e-1;

% Plot it with violins. Figure 6 in the manuscript

figure('position', [0, 0, 1150, 350], 'color', 'w')
[ha, pos] = tight_subplot(1,4,[.01 .01],[.175 .175],[.1 .1]); 
    % Third element is for width, and fourth is for height.
    
violincolor = [0.7, 0.7, 0.7];
markercolor = [0.7, 0.7, 0.7];

letters = {'a', 'b', 'c', 'd'}; 

for i=1:4
    axes(ha(i));
    
    % Note: the first dimension of prelambda is the site, the second is the
    % age group, and the third are the iterations.
    modeled = reshape(prelambda(:,i, :)*1e5, [size(prelambda(:,i, :),1), size(prelambda(:,i, :),3)])';
    
    modeled(:,isnan(logpointest(:,i)))=[];
    obs = logpointest(:,i);
    obs(isnan(logpointest(:,i))) = [];
    low = loglowerci(:,i);
    low(isnan(loglowerci(:,i))) = [];
    hi = logupperci(:,i);
    hi(isnan(logupperci(:,i))) = [];
    
    % store that to put it in Github
    % The code I put in github starts here
    errorbarxy(obs, median(log10(modeled))', ...
                low, hi, ...
                0,0, {'k.', [0.65 0.65 0.65], [0.65 0.65 0.65]})
    hold on
    violin(log10(modeled),'x', obs', 'mc', [], 'medc', [], 'bw', .25, ...
        'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', 1/2);  
    %set(gca, 'yscale', 'log','xscale', 'log')

    axis(log10([5e-1 1.5e4 5e-1 1.5e4]), 'square')

    if i<4
    jiggle = zeros(size(obs,1), 1);
    jiggle(obs==0) = (i<3)*0.075*(-1).^(1:sum((obs==0))); % so markers do not overlap with one anoether
    text(obs+jiggle, log10(median(modeled))', num2cell(find(~isnan(logpointest(:,i)))), 'horizontal','center', 'vertical','middle',  'BackgroundColor', [0.9 0.9 0.9], 'Margin', 1, 'fontWeight', 'Bold')
    else
    jiggle = 0.025*(-1).^find(~isnan(logpointest(:,i))); % so markers do not overlap with one anoether
    text(obs+jiggle, log10(median(modeled))', num2cell(find(~isnan(logpointest(:,i)))), 'horizontal','center', 'vertical','middle', 'BackgroundColor', [0.9 0.9 0.9], 'Margin', 1, 'fontWeight', 'Bold')        
    end
    
    refline(1,0)
    hold off
end

% text(log10(observed), log(median(modeled)'), num2cell(1:size(data_wide,1)), 'horizontal','left', 'vertical','middle')

minorticks = [log10(1:10), 1+log10(1:10), 2+log10(1:10), 3+log10(1:10), 4+log10(1:10)];

set(ha(1:4),'XTickLabel',''); set(ha(1:4),'YTickLabel','')
set(ha(1:4),'Xtick', 0:4); set(ha(1:4),'Ytick', 0:4)
set(ha(1:4),'XTickLabel',  {'1', '10', '100', '1,000', '10,000'}, 'XTickLabelRotation', 45, 'fontWeight', 'Bold')
set(ha(1),'YTickLabel', {'1', '10', '100', '1,000', '10,000'}, 'fontWeight', 'Bold')
% set(ha(1:5),'XMinorTick','on','YMinorTick','on')
% set(ha(1:5), 'XMinorTickValues', minorticks, 'YMinorTickValues', minorticks)

[ax, h1] = suplabel('Observed Cases per 100,000 Person-Years', 'x', [0.06 0.11 0.8800 0.8800]);
set(h1,'FontSize', 14, 'fontWeight', 'Bold')
[ax, h2a] = suplabel('Model Predicted Cases', 'y', [0.05 0.06 0.8800 0.8800]);
[ax, h2b] = suplabel('per 100,000 Person-Years', 'y', [0.065 0.06 0.8800 0.8800]);
set(h2a,'FontSize', 14, 'fontWeight', 'Bold')
set(h2b,'FontSize', 14, 'fontWeight', 'Bold')
[ax, h3] = suplabel('Model-Predicted vs Observed Cases', 't', [0.0100 0.025 0.9800 0.8800]);
set(h3,'FontSize', 18, 'fontWeight', 'Bold')

agegroups = {'Ages <2', 'Ages 2-4', 'Ages 5-14', 'Ages \geq 15'};

for i=1:4
title(ha(i), agegroups{1,i}, 'FontWeight', 'bold')
end
