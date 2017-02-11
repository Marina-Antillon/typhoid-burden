% Prepare the data as a long file
% prep_data_long

data.Age_Group = categorical(data.Age_Group);
% data.sel(strcmp(data.sel, 'KalkajiA 1996')) = {'Kalkaji 1996A'};
% data.sel(strcmp(data.sel, 'KalkajiB 1996')) = {'Kalkaji 1996B'};
data.sel = categorical(data.sel);
data.country = categorical(data.country);

data.survtype1 = categorical(data.survtype1);
% data.survtype2 = categorical(data.survtype2);

data.surveillance1 = zeros(size(data, 1), 1);
data.surveillance1(data.survtype1=='P') = 1; %  & data.survtype2=='P'
% data.surveillance2 = zeros(size(data, 1), 1);
% data.surveillance2(data.survtype1=='P' & data.survtype2=='H') = 1;

data.ages = nan(size(data,1), 1);
data.ages(data.Age_Group=='all') = 9;
data.ages(data.Age_Group=='Young Children 0-1') = 1;
data.ages(data.Age_Group=='Pre-K 2-4') = 2;
data.ages(data.Age_Group=='School 5-14') = 3;
data.ages(data.Age_Group=='Adults 15+') = 4;

data.ages(data.Age_Group=='Young & Pre-K 0-4') = 5;
data.ages(data.Age_Group=='Pre-K & School 2-14') = 6;
data.ages(data.Age_Group=='Young & Pre-K & School 0-14') = 7;
data.ages(data.Age_Group=='School & Adults 5+') = 8;

data.ages = categorical(data.ages);

data.cases = round(data.cases);
data.pop_time = round(data.pop_time);
data.offset = data.pop_time; 

data.GBD_super_regions_labels = data.GBD_super_regions;
data.GBD_super_regions = double(categorical(data.GBD_super_regions)); %% still need this

data.GBD_2015_labels = categorical(strrep(data.GBD_2015, '''', '')); %% still need this
data.GBD_2015 = double(categorical(strrep(data.GBD_2015, '''', ''))); %% still need this