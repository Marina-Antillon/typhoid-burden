% This is the routine that changes a long dataset to a wide dataset in
% order to simulate incidence and plot model-estimated incidence.

% Adjust the cases for the % people who didn't have a blood draw, etc!
% Assume NAN = 1.
data.part002(isnan(data.part002)) = 1;
data.part205(isnan(data.part205)) = 1;
data.part515(isnan(data.part515)) = 1;
data.part15o(isnan(data.part15o)) = 1;

data_wide = data(data.ages=='9',:); 

unique_sel = unique(data.sel);
data_wide.offseta = NaN(size(data_wide,1),1);
data_wide.offsetb = NaN(size(data_wide,1),1);
data_wide.offsetc = NaN(size(data_wide,1),1);
data_wide.offsetd = NaN(size(data_wide,1),1);
data_wide.offsete = NaN(size(data_wide,1),1);
data_wide.offsetf = NaN(size(data_wide,1),1);
data_wide.offsetg = NaN(size(data_wide,1),1);
data_wide.offseth = NaN(size(data_wide,1),1);
data_wide.offseti = NaN(size(data_wide,1),1);

data_wide.casesa = NaN(size(data_wide,1),1);
data_wide.casesb = NaN(size(data_wide,1),1);
data_wide.casesc = NaN(size(data_wide,1),1);
data_wide.casesd = NaN(size(data_wide,1),1);
data_wide.casese = NaN(size(data_wide,1),1);
data_wide.casesf = NaN(size(data_wide,1),1);
data_wide.casesg = NaN(size(data_wide,1),1);
data_wide.casesh = NaN(size(data_wide,1),1);
data_wide.casesi = NaN(size(data_wide,1),1);

data_wide.casesi = data_wide.cases;
data_wide.offseti = data_wide.offset;

for i=1:length(unique_sel)
    data_wide.offseta(data_wide.sel==unique_sel(i,1)) = data.pop_time(data.ages=='1' & data.sel==unique_sel(i,1)); 
    data_wide.offsetb(data_wide.sel==unique_sel(i,1)) = data.pop_time(data.ages=='2' & data.sel==unique_sel(i,1)); 
    data_wide.offsetc(data_wide.sel==unique_sel(i,1)) = data.pop_time(data.ages=='3' & data.sel==unique_sel(i,1));     
    data_wide.offsetd(data_wide.sel==unique_sel(i,1)) = data.pop_time(data.ages=='4' & data.sel==unique_sel(i,1));     
    data_wide.offsete(data_wide.sel==unique_sel(i,1)) = data.pop_time(data.ages=='5' & data.sel==unique_sel(i,1));     
    data_wide.offsetf(data_wide.sel==unique_sel(i,1)) = data.pop_time(data.ages=='6' & data.sel==unique_sel(i,1));     
    data_wide.offsetg(data_wide.sel==unique_sel(i,1)) = data.pop_time(data.ages=='7' & data.sel==unique_sel(i,1));     
    data_wide.offseth(data_wide.sel==unique_sel(i,1)) = data.pop_time(data.ages=='8' & data.sel==unique_sel(i,1));     
    
    data_wide.casesa(data_wide.sel==unique_sel(i,1)) = data.cases(data.ages=='1' & data.sel==unique_sel(i,1)); 
    data_wide.casesb(data_wide.sel==unique_sel(i,1)) = data.cases(data.ages=='2' & data.sel==unique_sel(i,1)); 
    data_wide.casesc(data_wide.sel==unique_sel(i,1)) = data.cases(data.ages=='3' & data.sel==unique_sel(i,1)); 
    data_wide.casesd(data_wide.sel==unique_sel(i,1)) = data.cases(data.ages=='4' & data.sel==unique_sel(i,1)); 
    data_wide.casese(data_wide.sel==unique_sel(i,1)) = data.cases(data.ages=='5' & data.sel==unique_sel(i,1)); 
    data_wide.casesf(data_wide.sel==unique_sel(i,1)) = data.cases(data.ages=='6' & data.sel==unique_sel(i,1)); 
    data_wide.casesg(data_wide.sel==unique_sel(i,1)) = data.cases(data.ages=='7' & data.sel==unique_sel(i,1)); 
    data_wide.casesh(data_wide.sel==unique_sel(i,1)) = data.cases(data.ages=='8' & data.sel==unique_sel(i,1)); 
end

for i = 1:size(data_wide,1)
if ~isnan(data_wide.casese(i))
    data_wide.offseta(i) = data_wide.popa(i)/(data_wide.popa(i) + data_wide.popb(i)).*data_wide.offsete(i);
    data_wide.offsetb(i) = data_wide.popb(i)/(data_wide.popa(i) + data_wide.popb(i)).*data_wide.offsete(i);
elseif ~isnan(data_wide.casesf(i))
    data_wide.offsetb(i) = data_wide.popb(i)/(data_wide.popb(i) + data_wide.popc(i)).*data_wide.offsetf(i);
    data_wide.offsetc(i) = data_wide.popc(i)/(data_wide.popb(i) + data_wide.popc(i)).*data_wide.offsetf(i);
elseif ~isnan(data_wide.casesg(i))
    data_wide.offseta(i) = data_wide.popa(i)/(data_wide.popa(i) + data_wide.popb(i) + data_wide.popc(i)).*data_wide.offsetg(i);
    data_wide.offsetb(i) = data_wide.popb(i)/(data_wide.popa(i) + data_wide.popb(i) + data_wide.popc(i)).*data_wide.offsetg(i);    
    data_wide.offsetc(i) = data_wide.popc(i)/(data_wide.popa(i) + data_wide.popb(i) + data_wide.popc(i)).*data_wide.offsetg(i);       
elseif ~isnan(data_wide.casesh(i))
    data_wide.offsetc(i) = data_wide.popc(i)/(data_wide.popc(i) + data_wide.popd(i)).*data_wide.offseth(i);
    data_wide.offsetd(i) = data_wide.popd(i)/(data_wide.popc(i) + data_wide.popd(i)).*data_wide.offseth(i);    
elseif ~isnan(data_wide.casesi(i))
    data_wide.offseta(i) = data_wide.popa(i)/100.*data_wide.offseti(i); 
    data_wide.offsetb(i) = data_wide.popb(i)/100.*data_wide.offseti(i); 
    data_wide.offsetc(i) = data_wide.popc(i)/100.*data_wide.offseti(i);
    data_wide.offsetd(i) = data_wide.popd(i)/100.*data_wide.offseti(i); 
end
end

data_wide.casese(or(~isnan(data_wide.casesa), ~isnan(data_wide.casesb))) = NaN; 
data_wide.casesf(or(~isnan(data_wide.casesb), ~isnan(data_wide.casesc))) = NaN; 
data_wide.casesg(or(or(~isnan(data_wide.casesa), ~isnan(data_wide.casesb)),~isnan(data_wide.casesc))) = NaN; 
data_wide.casesh(or(~isnan(data_wide.casesc), ~isnan(data_wide.casesd))) = NaN; 
data_wide.casesi(or(or(or(~isnan(data_wide.casesa), ~isnan(data_wide.casesb)), ~isnan(data_wide.casesc)), ~isnan(data_wide.casesd))) = NaN; 
