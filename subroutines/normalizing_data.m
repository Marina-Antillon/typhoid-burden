%% Normalizing
nrm = @(x) (x-nanmean(x(data.ages=='9')))/sqrt(nanvar(x(data.ages=='9')));
    % There are 9 copies of each value. You only want to normalize it once.
logit = @(x) (log(x./(100-x)));

data.wtr_ctr = nrm(logit(data.pipedwater));
data.san_ctr = nrm(logit(data.flushtoilets));
data.gini_ctr = nrm(data.gini_data); % the logit has a horrible outlier
data.ext_pov_ctr = nrm(logit(data.pov_hcr_200ppp_perc));
data.road_paved_ctr = nrm(logit(data.road_paved_perc)); % should it have been logit (TINY difference)
data.eduyrs_f_log_ctr = nrm(log(data.wed));
data.stunting_ctr = nrm(logit(data.stunting)); 
data.water_stress_log_ctr = nrm(log(data.water_blue));
data.popdens_log_ctr = nrm(log(data.popdens));
data.GDPcap_log_ctr = nrm(log(data.GDPcap));
data.floodrisk_log_ctr = nrm(log(data.floodrisk));
data.hivprev_ctr = nrm(log(data.hivprev)); 

% Record the normalization statistics.
% All this affects mapping and out-of-sample validations
% Must run it for LOO as well...

% If I change any variables above TRIPLE CHECK that these are correct 

normalize.wtr_mean = nanmean(logit(data.pipedwater(data.ages=='9')));
normalize.wtr_std = nanvar(logit(data.pipedwater(data.ages=='9')))^0.5;

normalize.san_mean = nanmean(logit(data.flushtoilets(data.ages=='9')));
normalize.san_std = nanvar(logit(data.flushtoilets(data.ages=='9')))^0.5;

normalize.gini_mean = nanmean(data.gini_data(data.ages=='9'));
normalize.gini_std = nanvar(data.gini_data(data.ages=='9'))^0.5;

normalize.ext_pov_mean = nanmean(logit(data.pov_hcr_200ppp_perc(data.ages=='9'))); 
normalize.ext_pov_std = nanvar(logit(data.pov_hcr_200ppp_perc(data.ages=='9')))^0.5;

normalize.road_paved_mean = nanmean(logit(data.road_paved_perc(data.ages=='9')));
normalize.road_paved_std = nanvar(logit(data.road_paved_perc(data.ages=='9')))^0.5;

normalize.eduyrs_f_log_mean = nanmean(log(data.wed(data.ages=='9')));
normalize.eduyrs_f_log_std = nanvar(log(data.wed(data.ages=='9')))^0.5;

normalize.stunting_mean = nanmean(logit(data.stunting(data.ages=='9')));
normalize.stunting_std = nanvar(logit(data.stunting(data.ages=='9')))^0.5;

normalize.GDPcap_log_mean = nanmean(log(data.GDPcap(data.ages=='9'))); 
normalize.GDPcap_log_std = nanvar(log(data.GDPcap(data.ages=='9')))^0.5; 

normalize.popdens_log_mean = nanmean(log(data.popdens(data.ages=='9')));
normalize.popdens_log_std = nanvar(log(data.popdens(data.ages=='9')))^0.5;

normalize.floodrisk_log_mean = nanmean(log(data.floodrisk(data.ages=='9')));
normalize.floodrisk_log_std = nanvar(log(data.floodrisk(data.ages=='9')))^0.5;

normalize.water_stress_log_mean = nanmean(log(data.water_blue(data.ages=='9')));
normalize.water_stress_log_std = nanvar(log(data.water_blue(data.ages=='9')))^0.5;

normalize.hivprev_mean = nanmean(log(data.hivprev(data.ages=='9')));
normalize.hivprev_std = nanvar(log(data.hivprev(data.ages=='9')))^0.5;

covar = {'popdens_log_ctr', 'GDPcap_log_ctr', 'gini_ctr', 'wtr_ctr', 'san_ctr', 'eduyrs_f_log_ctr', ...
         'road_paved_ctr', 'ext_pov_ctr', 'stunting_ctr', 'floodrisk_log_ctr', ...
         'hivprev_ctr', 'water_stress_log_ctr'};

% record max and min for each

normalize.wtr_min = min(data.pipedwater);
normalize.wtr_max = max(data.pipedwater); 

normalize.san_min = min(data.flushtoilets);
normalize.san_max = max(data.flushtoilets);

normalize.gini_min = min(data.gini_data);
normalize.gini_max = max(data.gini_data);

normalize.ext_pov_min = min(data.pov_hcr_200ppp_perc);
normalize.ext_pov_max = max(data.pov_hcr_200ppp_perc);

normalize.road_paved_min = min(data.road_paved_perc);
normalize.road_paved_max = max(data.road_paved_perc);

normalize.GDPcap_min = min(data.GDPcap); 
normalize.GDPcap_max = max(data.GDPcap); 

normalize.eduyrs_f_min = min(data.wed);
normalize.eduyrs_f_max = max(data.wed);

normalize.stunting_min = min(data.stunting);
normalize.stunting_max = max(data.stunting);

normalize.floodrisk_min = min(data.floodrisk);
normalize.floodrisk_max = max(data.floodrisk);

normalize.hivprev_min = min(data.hivprev);
normalize.hivprev_max = max(data.hivprev);

normalize.popdens_min = min(data.popdens);
normalize.popdens_max = max(data.popdens);

normalize.water_stress_min = min(data.water_blue);
normalize.water_stress_max = max(data.water_blue);
