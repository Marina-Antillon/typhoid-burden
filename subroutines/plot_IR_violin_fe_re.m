
% Condense the data for combinations of age groups 
combs = table2array(data_wide(:, {'casese', 'casesf', 'casesg', 'casesh', 'casesi'}));
combsoffset = table2array(data_wide(:, {'offsete', 'offsetf', 'offsetg', 'offseth', 'offseti'}));
[a, b] = find(~isnan(combs));

for j=1:2
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

figure('position', [0, 0, 1000, 450], 'color', 'w')
    % Third element is for width, and fourth is for height.
[ha, pos] = tight_subplot(2,5,[.01 .01],[.125 .125],[.075 .075]); 

for j=1:2
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

set(ha(1:10),'XTickLabel',''); set(ha(1:10),'YTickLabel','')
set(ha(1:10),'Xtick', 0:4); set(ha(1:10),'Ytick', 0:4)
set(ha(6:10),'XTickLabel',  {'1', '10', '100', '1,000', '10,000'}, 'XTickLabelRotation', 45)
set(ha([1,6]),'YTickLabel', {'1', '10', '100', '1,000', '10,000'})
% set(ha(1:15),'XMinorTick','on','YMinorTick','on')
% set(ha(1:15), 'XMinorTickValues', minorticks, 'YMinorTickValues', minorticks)
agegroups = {'Ages 0-2', 'Ages 2-5', 'Ages 5-15', 'Ages 15+', 'Combined Age Groups'};

for i=1:5
title(ha(i), agegroups{1,i}, 'FontWeight', 'normal')
end

[ax, h1] = suplabel('Observed cases per 100,000 people', 'x', [0.06 0.07 0.8800 0.8800]);
set(h1,'FontSize', 16)
[ax, h2a] = suplabel('Model-predicted cases per 100,000 people ', 'y', [0.05 0.06 0.8800 0.8800]);
% [ax, h2b] = suplabel('per 100,000 people', 'y', [0.045 0.06 0.8800 0.8800]);
set(h2a,'FontSize', 16)
% set(h2b,'FontSize', 24)
[ax, h3] = suplabel('Model-predicted vs observed cases', 't', [0.0500 0.07 0.9000 0.8800])
set(h3,'FontSize', 22)
