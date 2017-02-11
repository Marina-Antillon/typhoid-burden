combs = table2array(data_wide(:, {'casese', 'casesf', 'casesg', 'casesh', 'casesi'}));
combsoffset = table2array(data_wide(:, {'offsete', 'offsetf', 'offsetg', 'offseth', 'offseti'}));
[a, b] = find(~isnan(combs)); % find combination of age groups.

for i=1:length(a)
prelambdacombs(i, :) = prelambda(a(i), b(i)+4, :);
obscomb(i) = combs(a(i), b(i)); 
popcomb(i) = combsoffset(a(i), b(i)); 
end
% Violin plot

violincolor = [0.5, 0.5, 0.5];
markercolor = [0.75, 0 , 0];

letters = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'}; 

figure('position', [0, 0, 1500, 350], 'color', 'w')
[ha, pos] = tight_subplot(1,5,[.01 .01],[.1 .1],[.075 .075]); 
    % Third element is for width, and fourth is for height.

for i=1:5
    
if i<5
    axes(ha(i));
    observed = eval(['(data_wide.cases', letters{i}, './data_wide.offset', letters{i}, ')*1e5+1']);
    modeled = reshape(prelambda(:,i, :)*1e5, [size(prelambda(:,i, :),1), size(prelambda(:,i, :),3)])';

    bubblesize = eval(['max(1,data_wide.offset', letters{i}, '/1000).^0.75;']);
    
    modeled(:,isnan(observed))=[];
    bubblesize(isnan(observed))=[];
    observed(isnan(observed))=[]; 

    violin(log10(modeled),'x', log10(observed)', 'mc', [], 'medc', [], 'bw', 0.2, ...
    'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', .4);  
    %set(gca, 'yscale', 'log','xscale', 'log')
    hold on
    scatter(log10(observed), median(log10(modeled))', bubblesize, 'filled', ...
    'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)

    axis(log10([5e-1 1.5e4 5e-1 1.5e4]), 'square')
    refline(1,0)
    hold off
else
    axes(ha(5));
    observed = obscomb./popcomb*1e5+1;
    modeled = prelambdacombs'*1e5;

    bubblesize = max(1, popcomb/1000).^0.75;
    
    modeled(:,isnan(observed))=[];
    bubblesize(isnan(observed))=[];
    observed(isnan(observed))=[]; 
    
    violin(log10(modeled),'x', log10(observed)', 'mc', [], 'medc', [],  'bw', 0.2,...
    'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', .4);  
    %set(gca, 'yscale', 'log','xscale', 'log')
    hold on
    scatter(log10(observed), median(log10(modeled)), bubblesize, 'filled', ...
    'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)
    axis(log10([5e-1 1.5e4 5e-1 1.5e4]), 'square')
    refline(1,0)
    hold off
    
end

end

% numlabs = {'1', '100', '1,000', '10,000', '100,000'};  
% set(gca, 'YTick', single(1:5), 'XTick', single(1:5), ...
%     'XTickLabels', numlabs, 'YTickLabels', numlabs, ...
%     'XTickLabelRotation', 45);

minorticks = [log10(1:10), 1+log10(1:10), 2+log10(1:10), 3+log10(1:10), 4+log10(1:10)];

set(ha(1:5),'XTickLabel',''); set(ha(1:5),'YTickLabel','')
set(ha(1:5),'Xtick', 0:5); set(ha(1:5),'Ytick', 0:5)
set(ha(1:5),'XTickLabel',  {'1', '10', '100', '1,000', '10,000', '100,000'}, 'XTickLabelRotation', 45, 'fontWeight', 'Bold')
set(ha(1),'YTickLabel', {'1', '10', '100', '1,000', '10,000', '100,000'}, 'fontWeight', 'Bold')
% set(ha(1:5),'XMinorTick','on','YMinorTick','on')
% set(ha(1:5), 'XMinorTickValues', minorticks, 'YMinorTickValues', minorticks)

agegroups = {'Ages <2', 'Ages 2-4', 'Ages 5-14', 'Ages ?15', 'Combined Age Groups'};

for i=1:5
title(ha(i), agegroups{1,i}, 'FontWeight', 'bold')
end

[ax, h1] = suplabel('Observed Cases per 100,000 Person-Years', 'x', [0.06 0.125 0.8800 0.8800]);
set(h1,'FontSize', 14, 'fontWeight', 'Bold')
[ax, h2a] = suplabel('Model-Predicted Cases', 'y', [0.0325 0.06 0.8800 0.8800]);
[ax, h2b] = suplabel('per 100,000 Person-Years', 'y', [0.045 0.06 0.8800 0.8800]);
set(h2a,'FontSize', 14, 'fontWeight', 'Bold')
set(h2b,'FontSize', 14, 'fontWeight', 'Bold')
[ax, h3] = suplabel('Model-Predicted vs Observed Incidence', 't', [0.0100 0.025 0.9800 0.8800]); 
set(h3,'FontSize', 18, 'fontWeight', 'Bold')
