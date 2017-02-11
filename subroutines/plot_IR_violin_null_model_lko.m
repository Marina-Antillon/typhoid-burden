
% Condense the data for combinations of age groups 
combs = table2array(data_wide(:, {'casese', 'casesf', 'casesg', 'casesh', 'casesi'}));
combsoffset = table2array(data_wide(:, {'offsete', 'offsetf', 'offsetg', 'offseth', 'offseti'}));
[a, b] = find(~isnan(combs));

clear prelambdacombs obscomb popcomb
for j=1:3
for i=1:length(a)
prelambdacombs{j}(i, :) = prelambda{j}(a(i), b(i)+4, :);
obscomb{j}(i) = combs(a(i), b(i)); 
popcomb{j}(i) = combsoffset(a(i), b(i)); 
end
end

% Violin plot
violincolor = [0.5, 0.5, 0.5];
markercolor = [0.75 0 0.25];

letters = {'a', 'b', 'c', 'd'}; 

figure('position', [0, 0, 1000, 900], 'color', 'w')
[ha, pos] = tight_subplot(3,5,[.01 .01],[.1 .1],[.075 .075]); 
    % Third element is for width, and fourth is for height.

for j=1:3
    for i=1:5
        if i<5
            axes(ha((j-1)*5+i));
            observed = eval(['(data_wide.cases', letters{i}, './data_wide.offset', letters{i}, ')*1e5+1']);
            
        if j==1
            modeled = reshape(prelambda{j}(:,i, :)*1e5, [size(prelambda{j}(:,i, :),1), size(prelambda{j}(:,i, :),3)])';
        elseif j==2
            modeled = reshape(prelambda{j}(:,i, :)*1e5, [size(prelambda{j}(:,i, :),1), size(prelambda{j}(:,i, :),3)])';
        else
            modeled = reshape(prelambda{j}(:,i, :)*1e5, [size(prelambda{j}(:,i, :),1), size(prelambda{j}(:,i, :),3)])';
        end
        
            bubblesize = eval(['max(1,data_wide.offset', letters{i}, '/1000).^0.75;']);
    
            modeled(:,isnan(observed))=[];
            bubblesize(isnan(observed))=[];
            observed(isnan(observed))=[];
    
            violin(log10(modeled),'x', log10(observed)', 'mc', [], 'medc', [], ...
                'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', 1/2);  
            %set(gca, 'yscale', 'log','xscale', 'log')
            hold on
            scatter(log10(observed), median(log10(modeled)), bubblesize, 'filled', ...
            'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)

            axis(log10([5e-1 1.5e4 5e-1 1.5e4]), 'square')
            refline(1,0)
            hold off
        else
            axes(ha((j-1)*5+5));
            observed = obscomb{j}./popcomb{j}*1e5+1;
            modeled = prelambdacombs{j}'*1e5;

            bubblesize = max(1, popcomb{j}/1000).^0.75;
            
            modeled(:,isnan(observed))=[];
            bubblesize(isnan(observed))=[];
            observed(isnan(observed))=[]; 
            
            violin(log10(modeled),'x', log10(observed)', 'mc', [], 'medc', [], ...
            'facecolor', violincolor, 'edgecolor', [], 'facealpha', 0.3, 'scale', 1/2);  
            %set(gca, 'yscale', 'log','xscale', 'log')
            hold on
            scatter(log10(observed), median(log10(modeled)), bubblesize, 'filled', ...
            'MarkerFaceColor', markercolor, 'MarkerEdgeColor', markercolor)
            axis(log10([5e-1 1.5e4 5e-1 1.5e4]), 'square')
            refline(1,0)
            hold off    
        end
    end
end



minorticks = [log10(1:10), 1+log10(1:10), 2+log10(1:10), 3+log10(1:10), 4+log10(1:10)];

set(ha(1:15),'XTickLabel',''); set(ha(1:15),'YTickLabel','')
set(ha(1:15),'Xtick', 0:4); set(ha(1:15),'Ytick', 0:4)
set(ha(11:15),'XTickLabel',  {'1', '10', '100', '1,000', '10,000'}, 'XTickLabelRotation', 45, 'fontWeight', 'Bold')
set(ha([1,6,11]),'YTickLabel', {'1', '10', '100', '1,000', '10,000'}, 'fontWeight', 'Bold')
% set(ha(1:15),'XMinorTick','on','YMinorTick','on')
% set(ha(1:15), 'XMinorTickValues', minorticks, 'YMinorTickValues', minorticks)
agegroups = {'Ages <2', 'Ages 2-4', 'Ages 5-14', 'Ages ?15', 'Combined Age Groups'};

for i=1:5
title(ha(i), agegroups{1,i}, 'FontName','Arial', 'fontWeight', 'Bold')
end

[ax, h1] = suplabel('Observed Cases per 100,000 Person-Years', 'x', [0.06 0.05 0.8800 0.8800]);
set(h1,'FontSize', 16, 'fontWeight', 'Bold')

[ax, h2a] = suplabel('Model-Predicted Cases per 100,000 Person-Years', 'y', [0.05 0.06 0.8800 0.8800]);
% [ax, h2b] = suplabel('per 100,000 people', 'y', [0.045 0.06 0.8800 0.8800]);
set(h2a,'FontSize', 16, 'FontName','Arial', 'fontWeight', 'Bold')
% set(h2b,'FontSize', 24)
[ax, h3] = suplabel('Model-Predicted vs Observed Incidence', 't'); % , [0.0500 0.075 0.9000 0.8800]
set(h3,'FontSize', 22, 'FontName','Arial', 'fontWeight', 'Bold')

% text(-2, 0.9, plotletters{model}, 'FontSize', 14,'FontName','Arial') % , 'FontWeight', 'bold'

