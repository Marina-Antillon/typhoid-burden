% MAKING RASTER MAPS

% For the models where there is a regional component, I will bring in the
% dataset that partitions groups of countries:

% classifications = readtable('../MAP/world_map_regions_FINAL.txt', 'delimiter', '\t');
% 
% map_legend.WHO_region = getlevels(nominal(classifications.WHO_region));
% classifications.WHO_region = double(nominal(classifications.WHO_region));
% 
% map_legend.WHO_epi_region = getlevels(nominal(classifications.WHO_epi_region));
% classifications.WHO_epi_region = double(nominal(classifications.WHO_epi_region));
% 
% map_legend.GAVI_co_financing_DoV= getlevels(nominal(classifications.GAVI_co_financing_DoV));
% classifications.GAVI_co_financing_DoV = double(nominal(classifications.GAVI_co_financing_DoV));
% 
% map_legend.GBD_2015 = getlevels(nominal(classifications.GBD_2015));
% classifications.GBD_2015 = double(nominal(classifications.GBD_2015));
% 
% map_legend.GBD_super_regions = getlevels(nominal(classifications.GBD_super_regions));
% classifications.GBD_super_regions = double(nominal(classifications.GBD_super_regions));
% 
% map_legend.region_number = getlevels(nominal(classifications.region_number));
% classifications.region_number = double(nominal(classifications.region_number));

cd('/Users/Marina/Desktop/typhoid-burden/05-maps')

load('../00-data/world_map.mat') % the shp file, not a raster dataset.

load('../00-data/floods_highres.mat') 
load('../00-data/GDP_CAP_raster_all.mat') % TODO: update this to 2015
load('../00-data/mapped_wdi_vars.mat')

for j = 1:10
mapped_vars_wdi.(wdi_covar{j+2}) = mapped_vars{j,1};
end
% mapped_vars = mapped_vars_named;

load('../00-data/pop_mapres.mat')
load('../00-data/regions_indicator.mat')
load('../00-data/country_raster.mat')
load('../00-data/blue_water.mat')
load('../00-data/mapped_vars_gdl.mat')
load('../00-data/MASK2.mat')
load('../00-data/mapped_age_dist.mat')

% gdlmap = shaperead('../00-data/GDL_data/gdlplus.shp');

%% Logistic regression equation
logit = @(x) (log(x./(100-x)));

%% Population under study

for i=1:size(world,1)
population(i, 1) = world(i,1).pop_est;
population(i, 2) = world(i,1).GBD_2015;
end

population(population(:,1)==-99, :) = [];
population(ismember(population(:,2), [2, 5, 9, 11, 12, 20]), :) = [];

sum(population(:,1)) % 5.4947e+09

% Population in countries that have been LMIC at least once in the
% last 5 years.

geoshow(log(popfit), Rfinal, 'DisplayType','surface')
geoshow(log(gdp14cap+1), R_gdp, 'DisplayType','surface')

%% Predictors
% WDI & GDL vars - will need wdi_covar to pull out the right raster matrix
preddesc.san_ctr  = [nanmean(mapped_vars.flushtoilets(:).*mask(:)), prctile(mapped_vars.flushtoilets(:).*mask(:), [2.5, 100])];
preddesc.wtr_ctr = [nanmean(mapped_vars.pipedwater(:)), prctile(mapped_vars.pipedwater(:).*mask(:), [2.5, 100])];
preddesc.gini_ctr = [nanmean(mapped_vars_wdi.gini(:)), prctile(mapped_vars_wdi.gini(:).*mask(:), [2.5, 100])];
preddesc.ext_pov_ctr = [nanmean(mapped_vars_wdi.pov_hcr_200ppp_perc(:).*mask(:)), prctile(mapped_vars_wdi.pov_hcr_200ppp_perc(:).*mask(:), [2.5, 100])];
preddesc.eduyrs_f_log_ctr = [nanmean(log(mapped_vars.wed(:))), prctile(mapped_vars.wed(:).*mask(:), [2.5, 100])]; % LOG
preddesc.road_paved_log_ctr = [nanmean(mapped_vars_wdi.road_paved_perc(:).*mask(:)), prctile(mapped_vars_wdi.road_paved_perc(:).*mask(:), [2.5, 100])]; 

% Spatial dataset vars

min(min(popfit(popfit>0)))

preddesc.popdens_log_ctr = [exp(nanmean(log(popfit(:).*mask(:)+1))), prctile(popfit(:).*mask(:)+1, [2.5 100])]; % LOG
preddesc.GDPcap_log = [exp(nanmean(log(gdp14cap(:).*mask(:)+1))), prctile(gdp14cap(:).*mask(:), [2.5 100])]; % LOG 
preddesc.floodrisk_log_ctr = [exp(nanmean(log(floodsfit(:).*mask(:)+0.5))), prctile(floodsfit(:).*mask(:), [2.5 100])] ; % LOG
preddesc.water_stress_log_ctr = [exp(nanmean(log(bluemean(:).*mask(:)+0.5))), prctile(bluemean(:).*mask(:), [2.5 100])]; % LOG 

preddesc.stunting_ctr = [nanmean(mapped_vars.stunting(:).*mask(:)), prctile(mapped_vars.stunting(:).*mask(:), [2.5 100])];
preddesc.hivprev_ctr = [nanmean(mapped_vars.hivprev(:).*mask(:)), prctile(mapped_vars.hivprev(:).*mask(:), [2.5 100])]; %LOG

%                  san_ctr: [56.7720 0.3600 98.0034]
%                  wtr_ctr: [61.2326 1.2263 99.5386]
%                 gini_ctr: [26.2093 1.2263 99.5386]
%              ext_pov_ctr: [24.9802 0.0500 87.8700]
%         eduyrs_f_log_ctr: [1.7372 0.5975 10.6700]
%       road_paved_log_ctr: [46.8887 1.8200 91.1000]
%          popdens_log_ctr: [7.0634 1.0079 430.2667]
%               GDPcap_log: [6.9092e+03 694.6673 5.3836e+04]
%        floodrisk_log_ctr: [3.7221 0 32]
%     water_stress_log_ctr: [0.8461 3.0000e-04 434.7529]
%              hivprev_ctr: [1.4290 0.1000 13.3000]
%             stunting_ctr: [17.1446 0.0164 57.6000]

%% Normalizing the raster data according to the distribution parameters of the predictor data in the estimation sample.
% Lim to dataset min and max

load('../output/jan06_4ages_location_hp_slope.mat') 
% load('../spike_slab/normalize_fix.mat')

structArray(:,2)= [];

mapped_vars.flushtoilets(mapped_vars.flushtoilets<normalize.san_min) = normalize.san_min;
mapped_vars.flushtoilets(mapped_vars.flushtoilets>normalize.san_max) = normalize.san_max;
pred.san_ctr  = (logit(mapped_vars.flushtoilets) - normalize.san_mean)/normalize.san_std;

mapped_vars.pipedwater(mapped_vars.pipedwater<normalize.wtr_min) = normalize.wtr_min;
mapped_vars.pipedwater(mapped_vars.pipedwater>normalize.wtr_max) = normalize.wtr_max;
pred.wtr_ctr = (logit(mapped_vars.pipedwater) - normalize.wtr_mean)/normalize.wtr_std;

mapped_vars_wdi.gini(mapped_vars_wdi.gini<normalize.gini_min) = normalize.gini_min;
mapped_vars_wdi.gini(mapped_vars_wdi.gini>normalize.gini_max) = normalize.gini_max;
pred.gini_ctr = ((mapped_vars_wdi.gini) - normalize.gini_mean)/normalize.gini_std;

mapped_vars_wdi.pov_hcr_200ppp_perc(mapped_vars_wdi.pov_hcr_200ppp_perc<normalize.ext_pov_min) = normalize.ext_pov_min;
mapped_vars_wdi.pov_hcr_200ppp_perc(mapped_vars_wdi.pov_hcr_200ppp_perc>normalize.ext_pov_max) = normalize.ext_pov_max;
pred.ext_pov_ctr = (logit(mapped_vars_wdi.pov_hcr_200ppp_perc.*mask) - normalize.ext_pov_mean)/normalize.ext_pov_std;

mapped_vars.wed(mapped_vars.wed<normalize.eduyrs_f_min) = normalize.eduyrs_f_min;
mapped_vars.wed(mapped_vars.wed>normalize.eduyrs_f_max) = normalize.eduyrs_f_max;
pred.eduyrs_f_log_ctr = (log(mapped_vars.wed) - normalize.eduyrs_f_log_mean)/normalize.eduyrs_f_log_std;

mapped_vars_wdi.road_paved_perc(mapped_vars_wdi.road_paved_perc<normalize.road_paved_min) = normalize.road_paved_min;
mapped_vars_wdi.road_paved_perc(mapped_vars_wdi.road_paved_perc>normalize.road_paved_max) = normalize.road_paved_max;
pred.road_paved_ctr = (logit(mapped_vars_wdi.road_paved_perc).*mask - normalize.road_paved_mean)/normalize.road_paved_std;

% stunting
mapped_vars.stunting(mapped_vars.stunting<normalize.stunting_min) = normalize.stunting_min;
mapped_vars.stunting(mapped_vars.stunting>normalize.stunting_max) = normalize.stunting_max;
pred.stunting_ctr = (logit(mapped_vars.stunting) - normalize.stunting_mean)/normalize.stunting_std;
% HIV
mapped_vars.hivprev(mapped_vars.hivprev<normalize.hivprev_min) = normalize.hivprev_min;
mapped_vars.hivprev(mapped_vars.hivprev>normalize.hivprev_max) = normalize.hivprev_max;
pred.hivprev_ctr = (log(mapped_vars.hivprev) - normalize.hivprev_mean)/normalize.hivprev_std;

% Spatial dataset vars
popfit2=popfit;

popfit2(popfit2<normalize.popdens_min) = normalize.popdens_min;
popfit2(popfit2>normalize.popdens_max) = normalize.popdens_max;
pred.popdens_log_ctr = (log(popfit2.*mask) - normalize.popdens_log_mean)/normalize.popdens_log_std;

gdp14cap(gdp14cap<normalize.GDPcap_min) = normalize.GDPcap_min; % HERE
gdp14cap(gdp14cap>normalize.GDPcap_max) = normalize.GDPcap_max;
pred.GDPcap_log_ctr = (log(gdp14cap.*mask+1) - normalize.GDPcap_log_mean)/normalize.GDPcap_log_std;

floodsfit(floodsfit<normalize.floodrisk_min) = normalize.floodrisk_min;
floodsfit(floodsfit>normalize.floodrisk_max) = normalize.floodrisk_max;
pred.floodrisk_log_ctr = (log(floodsfit.*mask+0.5) - normalize.floodrisk_log_mean)/normalize.floodrisk_log_std;

% water stress
bluemean(bluemean<normalize.water_stress_min) = normalize.water_stress_min;
bluemean(bluemean>normalize.water_stress_max) = normalize.water_stress_max;
pred.water_stress_log_ctr = (log(bluemean.*mask) - normalize.water_stress_log_mean)/normalize.water_stress_log_std;

%% Mask

oldmask = region_map{1,1}; 
oldmask(ismember(oldmask, [2, 5, 9, 11, 12, 20])) = nan;
oldmask(~isnan(oldmask)) = 1;

newmask = mask;
newmask(~isnan(oldmask) & isnan(newmask)) = 1;
geoshow(newmask, Rfinal, 'DisplayType','surface')

% % Area in each 1/10 degree squared, haversine formula
r = 6371;
lat = [90:-0.1:-89.9]'; 
% a = (sin(0.5*0.1*pi/180))^2 + cos(lat*pi/180).*cos((lat-0.1)*pi/180)*(sin(0.5*0.1*pi/180))^2;
% % c = 2*atan2(sqrt(a), sqrt(1-a));
% d = 2*r*asin(sqrt(a));

% area = 0.5*(repmat(d(1:end-1)', 1, 3599)+repmat(d(2:end)', 1, 3599))*pi*r*0.1/180; 
area_sphere = 4*pi*r^2;

% Theorem of Archimedes for a spherical model of the Earth

area = repmat((sin(lat(1:end-1)*pi/180) - sin(lat(2:end)*pi/180)) * (0.1*pi/180) * r^2, 1, 3600);
area = [area; area(end,:)];

all = nansum(nansum(popfit.*mask.*area)); % 6.04*10^9 people

all_old = nansum(nansum(popfit.*oldmask.*area)); % 5.8525e+09 or 6.51% more people than I had expected.

%% RASTER MAP

mask_vector = mask(:);
index_mask = find(~isnan(mask_vector));

for i = 1:length(covar)
    input.temp(:,i) = (pred.(covar{i})(:));
    % input.cov(:,i) = input.cov(:,i).*mask(:);
end

input.cov = input.temp(index_mask, :);
input = rmfield(input, 'temp');

% Reduce input.cov further to only the non-nan (non-nan according to mask) parts of the world.

clearvars -except input mask mask_vector index_mask structArray popfit popfit2 R_pop oldmask mapped_age_dist

tmp  =  single(mapped_age_dist{1,1}*1/6*1/100); 
age1dist = tmp(:);
age1dist = age1dist(index_mask);
tmp =  single(mapped_age_dist{1,1}*1/6*1/100); 
age2dist = tmp(:);
age2dist = age2dist(index_mask);
tmp =  single(mapped_age_dist{1,1}*2/3*1/100); 
age3dist = tmp(:);
age3dist = age3dist(index_mask);
tmp = single(mapped_age_dist{2,1}*1/100+mapped_age_dist{3,1}*1/100); 
age4dist = tmp(:);
age4dist = age4dist(index_mask);

cd('../output')
draws = 100;
total_rounds = round(1000/draws);

% save all in csv format...

% I split up the random draws because the computer couldn't take it
% altogether.

tic
for round = 1:total_rounds
round % show me what round this is in
post_samp = randsample(1:length(structArray(1,1).mu),draws);

% Age group 1 

pred_mat_all1 = single(nan(size(index_mask, 1), draws)); 

for i=1:draws

j = post_samp(i);
par.mu = (structArray(1,1).mu(j,:));
par.beta = (structArray(1,1).beta(j,:));
par.gamma = (structArray(1,1).gamma(j,:)); 
par.gamma_int1 = (structArray(1,1).gamma_int1(j,:));  

lambda = exp(par.mu + par.beta(1) + par.gamma*input.cov' + par.gamma_int1*input.cov' + log(1e5));

lambda(lambda<1) = 1;
lambda(lambda>1e5) = 1e5;

% pred_mat_all1(:,:,i) = single(reshape(lambda, size(popfit_exp, 1), size(popfit_exp, 2))); 

pred_mat_all1(:, i) = single(lambda');

end
pred_mat_all1 = single(pred_mat_all1);

filename = ['../output/predicted_map_age1_round', sprintf('%02.0f', round),'.csv'];

dlmwrite(filename, pred_mat_all1)

clear par pred_mat_all1

% Age group 2

pred_mat_all2 = single(nan(size(index_mask, 1), draws)); 

for i=1:draws

j = post_samp(i);
par.mu = (structArray(1,1).mu(j,:));
par.beta = (structArray(1,1).beta(j,:));
par.gamma = (structArray(1,1).gamma(j,:)); 
par.gamma_int2 = (structArray(1,1).gamma_int2(j,:));  

lambda = exp(par.mu + par.beta(2) + par.gamma*input.cov' + par.gamma_int2*input.cov' + log(1e5));

lambda(lambda<1) = 1;
lambda(lambda>1e5) = 1e5;

pred_mat_all2(:, i) = single(lambda');

end
pred_mat_all2 = single(pred_mat_all2);

filename = ['../output/predicted_map_age2_round', sprintf('%02.0f', round),'.csv'];

dlmwrite(filename, pred_mat_all2)

clear par pred_mat_all2

pred_mat_all3 = single(nan(size(index_mask, 1), draws)); 

for i=1:draws

j = post_samp(i);
par.mu = (structArray(1,1).mu(j,:));
par.gamma = (structArray(1,1).gamma(j,:)); 

lambda = exp(par.mu + par.gamma*input.cov' + log(1e5));

lambda(lambda<1) = 1;
lambda(lambda>1e5) = 1e5;

pred_mat_all3(:, i) = single(lambda');

end
pred_mat_all3 = single(pred_mat_all3);

filename = ['../output/predicted_map_age3_round', sprintf('%02.0f', round),'.csv'];

dlmwrite(filename, pred_mat_all3)

clear par pred_mat_all3

pred_mat_all4 = single(nan(size(index_mask, 1), draws)); 

for i=1:draws

j = post_samp(i);
par.mu = (structArray(1,1).mu(j,:));
par.beta = (structArray(1,1).beta(j,:));
par.gamma = (structArray(1,1).gamma(j,:)); 
par.gamma_int3 = (structArray(1,1).gamma_int3(j,:));  

lambda = exp(par.mu + par.beta(3) + par.gamma*input.cov' + par.gamma_int3*input.cov' + log(1e5));

lambda(lambda<1) = 1;
lambda(lambda>1e5) = 1e5;

pred_mat_all4(:, i) = single(lambda');
end
pred_mat_all4 = single(pred_mat_all4);

filename = ['../output/predicted_map_age4_round', sprintf('%02.0f', round),'.csv'];
dlmwrite(filename, pred_mat_all4)

clear par pred_mat_all4;

end
toc
% This took 2700 seconds

% Make the file for the overall incidence distribution
tic
for round=1:total_rounds
    for age = 1:4
    filename = ['../output/predicted_map_age', sprintf('%1.0f', age) '_round', sprintf('%02.0f', round),'.csv'];
    arrayname = ['pred_mat_all', sprintf('%1.0f', age)]; 
    eval([arrayname, ' = dlmread(filename);']);
    end

% Add the cases by age group.
pred_mat_all5 = pred_mat_all1.*repmat(age1dist, [1,size(pred_mat_all1,2)])...
    + pred_mat_all2.*repmat(age2dist, [1,size(pred_mat_all2,2)]) ...
    + pred_mat_all3.*repmat(age3dist, [1,size(pred_mat_all3,2)]) ...
    + pred_mat_all4.*repmat(age4dist, [1,size(pred_mat_all4,2)]);

pred_mat_all5 = single(pred_mat_all5);

filename = ['../output/predicted_map_age5_round', sprintf('%02.0f', round),'.csv'];
dlmwrite(filename, pred_mat_all5)

clear par pred_mat_all;

end
toc % 1665 seconds

% Concatenate the maps for each age group together

tic
for age = 1:5
    for round=1:total_rounds
    filename = ['../output/predicted_map_age', sprintf('%1.0f', age) '_round', sprintf('%02.0f', round),'.csv'];
    objname_new{round, 1} = dlmread(filename);
    end

% Concatenate the rounds
arrayname = ['pred_mat_age', sprintf('%1.0f', age)];
if total_rounds>1
    eval([arrayname, ' = cat(2, objname_new{1,1}, objname_new{2,1});']);
else
     eval([arrayname, ' = objname_new{1,1} ;'])
end

if total_rounds>2
    for round=3:total_rounds
    eval([arrayname, ' = cat(2,', arrayname ', objname_new{round,1} );']);
    end
end

filename = ['../output/predicted_map_age', sprintf('%1.0f', age), '_concatenated.csv'];
eval(['dlmwrite(filename, ' arrayname, ');'])

eval(['clear ', arrayname, ' objname_new']);
end
toc % 5358 seconds

clear pred_mat_all1 pred_mat_all2 pred_mat_all3 pred_mat_all1 pred_mat_all4 pred_mat_all5;

% SUMMARY STATISTICS

% Median for each age group. 
tic
for age=1:5
filename = ['../output/predicted_map_age', sprintf('%1.0f', age), '_concatenated.csv'];
arrayname = ['pred_mat_age', sprintf('%1.0f', age)];

eval([arrayname, ' = dlmread(filename);']);
    
pred_medians(:,age) = eval(['median(', arrayname, ',2)']); 
eval(['clear ', arrayname]);
end
toc % 1728

pred_medians(pred_medians<1) = 1;
pred_medians(pred_medians>10000) = 9999;
pred_medians = single(pred_medians);
save('../output/pred_medians.mat', 'pred_medians', 'R_pop', '-v7.3')

% Probability of belonging to each incidence group.
pred_mat_age5 = dlmread('../output/predicted_map_age5_concatenated.csv');

pred_cat1 = mean((pred_mat_age5<10), 2);
pred_cat1(isnan(pred_cat1)) = NaN;
pred_cat2 = mean((pred_mat_age5>10 & pred_mat_age5<=100), 2); 
pred_cat2(isnan(pred_cat2)) = NaN;
pred_cat3 = mean((pred_mat_age5>100 & pred_mat_age5<=500), 2); 
pred_cat3(isnan(pred_cat3)) = NaN;
pred_cat4 = mean((pred_mat_age5>500), 2); 
pred_cat4(isnan(pred_cat4)) = NaN;

save('../output/pred_categories.mat', 'pred_cat1', 'pred_cat2', 'pred_cat3', 'pred_cat4', 'R_pop', '-v7.3')

%% Make a pretty map!

% load('./output/index_mask.mat') % WRONG. New one is nec. 
load('../00-data/pop_mapres.mat')

land = shaperead('landareas', 'UseGeoCoords', true);
% load('coast') brings in lat and long

load('../output/pred_medians.mat')
pred_medians_big = nan(1800*3600, 5);
pred_medians_big(index_mask,:) = pred_medians();

pred_med_map1 = single(reshape(pred_medians_big(:,1), size(popfit, 1), size(popfit, 2)));
pred_med_map2 = single(reshape(pred_medians_big(:,2), size(popfit, 1), size(popfit, 2)));
pred_med_map3 = single(reshape(pred_medians_big(:,3), size(popfit, 1), size(popfit, 2)));
pred_med_map4 = single(reshape(pred_medians_big(:,4), size(popfit, 1), size(popfit, 2)));
pred_med_map5 = single(reshape(pred_medians_big(:,5), size(popfit, 1), size(popfit, 2)));

% Overall incidence
figure
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,log10(10000)])
colormap((parula(100))) %flipud

ax = worldmap('World');
% setm(ax, 'Origin', [180 0 180])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(log10(pred_med_map5), R_pop, 'DisplayType','surface')

% plotm(lat, long)
title('Total Incidence per 100K persons', 'FontSize', 14)
c = {'1', '10', '100', '1,000', '10,000'};
colorbar('FontSize',11,'YTick',0:4,'YTickLabel',c);
[cb1, cb2] = caxis; % 1.9, 3.98; we are almost there...

% By age group:

colormap(parula(100))
figure('position', [0, 0, 1200, 500], 'color', 'w')
[ha, pos] = tight_subplot(2,2,[.01 .01],[.05 .05],[.15 .15]); 

axes(ha(1))
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,4])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(log10(pred_med_map1), R_pop, 'DisplayType','surface')
title({'Incidence per 100K persons'; 'ages 0-1'})
% colorbar
axes(ha(2))
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,4])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(log10(pred_med_map2), R_pop, 'DisplayType','surface')
title({'Incidence per 100K persons'; 'ages 2-4'})
% colorbar
axes(ha(3))
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,4])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(log10(pred_med_map3), R_pop, 'DisplayType','surface')
title({'Incidence per 100K persons'; 'ages 5-15'})
% colorbar
axes(ha(4))
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,4])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(log10(pred_med_map4), R_pop, 'DisplayType','surface')
title({'Incidence per 100K persons'; 'ages 15+'})

ha4p = ha(4).Position; % given as [x y width height], same as below
c = {'1', '10', '100', '1,000', '10,000'};
colorbar('FontSize',11,'YTick',0:4,'YTickLabel',c, 'Position', [ha4p(1)+ha4p(3)+0.025  ha4p(2)+0.1*ha4p(4)  0.01  ha4p(4)*1.8])

% Probability of belonging to each category

load('../output/pred_categories.mat')

pred_cat_big = nan(1800*3600, 4);
pred_cat_big(index_mask,1) = pred_cat1;
pred_cat_big(index_mask,2) = pred_cat2;
pred_cat_big(index_mask,3) = pred_cat3;
pred_cat_big(index_mask,4) = pred_cat4;

pred_cat_map1 = single(reshape(pred_cat_big(:,1), size(popfit, 1), size(popfit, 2)));
pred_cat_map2 = single(reshape(pred_cat_big(:,2), size(popfit, 1), size(popfit, 2)));
pred_cat_map3 = single(reshape(pred_cat_big(:,3), size(popfit, 1), size(popfit, 2)));
pred_cat_map4 = single(reshape(pred_cat_big(:,4), size(popfit, 1), size(popfit, 2)));

figure('position', [0, 0, 1200, 500], 'color', 'w')
[ha, pos] = tight_subplot(2,2,[.01 .01],[.05 .05],[.15 .15]); 

colormap(parula(100))
axes(ha(1))
caxis([0,1])
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(pred_cat_map1, R_pop, 'DisplayType','surface')
title({'Probability of Low Incidence:'; '<10 per 100,000 Person-Years'}) 

axes(ha(2))
caxis([0,1])
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(pred_cat_map2, R_pop, 'DisplayType','surface')
title({'Probability of Medium Incidence:'; '10-<100 per 100,000 Person-Years'})

axes(ha(3))
caxis([0,1])
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(pred_cat_map3, R_pop, 'DisplayType','surface')
title({'Probability of High Incidence:'; '100-<500 per 100,000 Person-Years'}) 

axes(ha(4))
caxis([0,1])
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
geoshow(pred_cat_map4, R_pop, 'DisplayType','surface')
title({'Probability of Very High Incidence:'; '500+ per 100,000 Person-Years'}) 

ha4p = ha(4).Position; % given as [x y width height], same as below
colorbar('Position', [ha4p(1)+ha4p(3)+0.025  ha4p(2)+0.1*ha4p(4)  0.01  ha4p(4)*1.8])

%% Calculating incidence.

load('../00-data/MASK2.mat')
mask_vector = mask(:);
index_mask = find(~isnan(mask_vector));

load('../MAP/pop_mapres.mat')
% pred_mat_age5 = dlmread('C:/Users/Marina Antillon/Google Drive/map_output/dec_map/predicted_map_age5_concatenated.csv');
pred_mat_age5 = dlmread('c:/Users/Marina Antillon/Google Drive/map_output/map_/predicted_map_age5_concatenated.csv'); 

load('../MAP/world_map.mat') % the shp file, not a raster dataset.
for i=1:size(world,1)
iso3{i,1} = world(i,1).ISO3;
% LAT(i,1) = world(i,1).Lat;
% LON(i,1) = world(i,1).Lon;
% REGION(i,1)=world(i,1).REGION;
end

% WORLDWIDE

% Theorem of Archimedes for a spherical model of the Earth
r = 6371;
lat = [90:-0.1:-89.9]'; 
area = repmat((sin(lat(1:end-1)*pi/180) - sin(lat(2:end)*pi/180)) * (0.1*pi/180) * r^2, 1, 3600);
area = [area; area(end,:)];

% index_mask_old = find(~isnan(oldmask(:)));

% index_mask_new2old = ismember(index_mask_old, index_mask);

pred_mat_age5(pred_mat_age5>10000) = 10000;
% pred_mat_age5(pred_mat_age5<1) = 1;

pop_redux = popfit(:);
pop_redux = pop_redux(index_mask);
area_redux = area(:);
area_redux = area_redux(index_mask);

% sum of the median (WRONG)
% all = nansum(nansum(pop_redux.*area_redux.*median(pred_mat_age5, 2)/(10^5))); 

% Median of the sum (CORRECT)
for i=1:size(pred_mat_age5, 2)
    cases_all(i) = nansum(pred_mat_age5(:,i).*pop_redux(:).*area_redux(:)/(10^5));
end

    all = nansum(pop_redux(:).*area_redux(:))

% for i=1:size(pred_mat_age5, 2)
%     cases_all(i) = nansum(pred_mat_age5(index_mask_new2old,i).*pop_redux(index_mask_new2old).*area_redux(index_mask_new2old)/(10^5));
% end
hist(cases_all/1e6,100)
prctile(cases_all, [2.5 50 97.5])./1e6 
prctile(cases_all(1:500), [2.5 50 97.5])./1e6
prctile(cases_all(501:1000), [2.5 50 97.5])./1e6 

% right before Gates: 8.1429   38.0408  130.7384

prctile(cases_all, [2.5 50 97.5])./1e6;
prctile(cases_all, [2.5 50 97.5])./all*1e5

save('../output/cases_all_jan06.mat', 'cases_all', 'all') 

% REGION ESTIMATES 

load('../MAP/regions_indicator.mat')

load('../MAP/GDL_data/MASK2.mat')
mask_vector = mask(:);
index_mask = find(~isnan(mask_vector));

region_redux = region_map{2,1}(:);
geoshow(region_map{2,1}, R_pop, 'DisplayType','surface')
region_redux = region_redux(index_mask);

region_study = nan(1800*3600, 1);
region_study(index_mask) = region_redux;
region_study_map = single(reshape(region_study, 1800, 3600));
geoshow(region_study_map, R_pop, 'DisplayType','surface')
hist(region_study_map(:), 7)

cases_gbdsr = nan(7, size(pred_mat_age5, 2));
incidence_gbdsr = nan(7, size(pred_mat_age5, 2));

for j=1:7
    tmp_pop = pop_redux(region_redux==j);
    tmp_area = area_redux(region_redux==j);
    tmp_pop_weight = (tmp_pop.*tmp_area)/nansum(tmp_pop.*tmp_area);
    for i=1:size(pred_mat_age5, 2)
        cases_gbdsr(j,i)=nansum(pred_mat_age5(region_redux==j,i).*tmp_pop.*tmp_area/(10^5));
        incidence_gbdsr(j,i)=nansum(pred_mat_age5(region_redux==j,i).*(10^5));
    end
end

region_labels = map_legend.GBD_super_regions';
save('../output/cases_inc_gbdsr_jan06.mat', 'cases_gbdsr', 'incidence_gbdsr', 'region_labels') 

cases_gbdSR_TABLE = cell2table(cellstr(region_labels), 'VariableNames', {'Region'});
cases_gbdSR_TABLE = [cases_gbdSR_TABLE, array2table(prctile(cases_gbdsr, [50 2.5 97.5], 2))]; 

writetable(cases_gbdSR_TABLE, '../output/cases_gbdSR_TABLE_jan06.csv')

incidence_gbdSR_TABLE = [cell2table(cellstr(region_labels), 'VariableNames', {'Region'}), array2table(prctile(incidence_gbdsr, [50 2.5 97.5], 2))]; 
writetable(incidence_gbdSR_TABLE, '../output/incidence_gbdSR_TABLE_jan06.csv')

% GBD REGIONS
% use region_map{3,1}

region_redux = region_map{1,1}(:);
geoshow(region_map{3,1}, R_pop, 'DisplayType','surface')
region_redux = region_redux(index_mask);

region_study = nan(1800*3600, 1);
region_study(index_mask) = region_redux;
region_study_map = single(reshape(region_study, 1800, 3600));
geoshow(region_study_map, R_pop, 'DisplayType','surface')
hist(region_study_map(:), 21)

cases_gbdreg = nan(21, size(pred_mat_age5, 2));
incidence_gbdreg = nan(21, size(pred_mat_age5, 2));

for j=1:21
    tmp_pop = pop_redux(region_redux==j);
    tmp_area = area_redux(region_redux==j);
    tmp_pop_weight = (tmp_pop.*tmp_area)/nansum(tmp_pop.*tmp_area);
    for i=1:size(pred_mat_age5, 2)
        cases_gbdreg(j,i)=nansum(pred_mat_age5(region_redux==j,i).*tmp_pop.*tmp_area/(10^5));
        incidence_gbdreg(j,i)=nansum(pred_mat_age5(region_redux==j,i).*tmp_pop_weight/(10^5));
    end
end

gbdregion_labels = map_legend.GBD_2015'; % labels
save('../output/cases_inc_gbdreg_jan06.mat', 'cases_gbdreg', 'incidence_gbdreg', 'gbdregion_labels') 

cases_gbdreg_TABLE = cell2table(cellstr(gbdregion_labels), 'VariableNames', {'Region'});
cases_gbdreg_TABLE = [cases_gbdreg_TABLE, array2table(prctile(cases_gbdreg, [50 2.5 97.5], 2))]; 
writetable(cases_gbdreg_TABLE, '../output/cases_gbdreg_TABLE_jan06.csv')

incidence_gbdreg_TABLE = [cell2table(cellstr(gbdregion_labels), 'VariableNames', {'Region'}), array2table(prctile(incidence_gbdreg, [50 2.5 97.5], 2))]; 
writetable(incidence_gbdreg_TABLE, '../output/incidence_gbdreg_TABLE_jan06.csv')

%% Probability that each region is in each of the incidence brackets.

% GBD region
cat_reg(:,1) = mean(single(incidence_gbdreg<10/1e5), 2)*100; 
cat_reg(:,2) = mean(single(10/1e5<=incidence_gbdreg & incidence_gbdreg<100/1e5), 2)*100; 
cat_reg(:,3) = mean(single(100/1e5<=incidence_gbdreg & incidence_gbdreg<500/1e5), 2)*100; 
cat_reg(:,4) = mean(single(500/1e5<=incidence_gbdreg), 2)*100; 

cat_reg_labels = cell2table(cellstr(gbdregion_labels), 'VariableNames', {'Region'});
cat_reg = [cat_reg_labels, array2table(cat_reg)];

% GBD region
cat_SR(:,1) = mean(single(incidence_gbdsr<10/1e5), 2)*100; 
cat_SR(:,2) = mean(single(10/1e5<=incidence_gbdsr & incidence_gbdsr<100/1e5), 2)*100; 
cat_SR(:,3) = mean(single(100/1e5<=incidence_gbdsr & incidence_gbdsr<500/1e5), 2)*100; 
cat_SR(:,4) = mean(single(500/1e5<=incidence_gbdsr), 2)*100;

cat_SR_labels = cell2table(cellstr(region_labels), 'VariableNames', {'Region'});
cat_SR = [cat_SR_labels, array2table(cat_SR)];

writetable(cat_SR, '../output/cat_SR_jan06.csv')
writetable(cat_reg, '../output/cat_reg_jan06.csv')

violincolor = [0.7, 0.7, 0.7];
markercolor = [0.7, 0.7, 0.7];

violin(([log10(cases_all'), log10(cases_gbd_superregion([1, 3:7],:)')]), ...
    'facecolor', violincolor, 'bw', 0.1, 'edgecolor', [], 'facealpha', 1,...
    'scale', 1/2, 'mc', [], 'medc', []);


set(gca, 'XTick', 1:7, 'XTickLabel', region_labels, 'XTickLabelRotation', 30,...
    'YTick', 4:8, 'YTickLabel', {'0.01', '0.1', '1', '10', '100'}, 'FontSize', 12)

%     GBD_super_regions                 GBD_super_regions_labels             
%     _________________    __________________________________________________
% 
%     1                    'Central Europe, Eastern Europe, and Central Asia'
%     2                    'Latin America and Caribbean'                     
%     3                    'North Africa and Middle East'                    
%     4                    'South Asia'                                      
%     5                    'Southeast Asia, East Asia, Oceania'              
%     6                    'Sub-Saharan Africa' 
    
    
% According to GBD 2015 classifications...
% {'Andean Latin America', 'Caribbean', 'Central Europe', 'Central Asia', 'Central Latin America', ...
% 'Central Sub-Saharan Africa', 'East Asia', 'Eastern Europe', 'Eastern Sub-Saharan Africa', ...
% 'North Africa and Middle East', 'South Asia', 'Southeast Asia', 'Southern Latin America', ...
% 'Southern Sub-Saharan Africa', 'Tropical Latin America', 'Western Sub-Saharan Africa'} 
% 16 but there are 14 in my map region_map from which mask is made... 
% Where is Oceania? Just complete out of the map?

