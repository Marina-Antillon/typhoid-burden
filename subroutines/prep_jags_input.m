% Preparing the input values for a jags run. 
% Works for either the GLM model or the markov model.

% take out overall incidence and incidence in <5 when it is not needed.
datajags=data;
datajags(or(isnan(datajags.cases), and(datajags.ages=='9', datajags.augmentation~=1)),:) = []; 
    % cases for the whole population are in the dataset, but this usually
    % repeats cases that are listed in records for age-specific incidence
datajags(or(isnan(datajags.cases), and(datajags.ages=='5', datajags.augmentation~=2)),:) = []; 
    % cases for the population <5 are in the dataset, but this usually
    % repeats cases that are listed in records for groups <2 and 2-4
    
% Adjust the cases for the % people who didn't have a blood draw, etc!
% Assume NAN = 1.
datajags.part002(isnan(datajags.part002)) = 1;
datajags.part205(isnan(datajags.part205)) = 1;
datajags.part515(isnan(datajags.part515)) = 1;
datajags.part15o(isnan(datajags.part15o)) = 1;

input.part(datajags.ages=='1') = datajags.part002(datajags.ages=='1');
input.part(datajags.ages=='2') = datajags.part205(datajags.ages=='2');
input.part(datajags.ages=='3') = datajags.part515(datajags.ages=='3');
input.part(datajags.ages=='4') = datajags.part15o(datajags.ages=='4');
input.part(datajags.ages=='5') = datajags.part205(datajags.ages=='5');
input.part(datajags.ages=='6') = datajags.part515(datajags.ages=='6');
input.part(datajags.ages=='7') = datajags.part515(datajags.ages=='7');
input.part(datajags.ages=='8') = datajags.part515(datajags.ages=='8');
input.part(datajags.ages=='9') = datajags.part515(datajags.ages=='9');

input.cases = datajags.cases';
input.ages = double(datajags.ages)';
input.surv = double(datajags.survtype1=='P')'; 
input.cov = table2array(datajags(:,covar));
input.country = double(categorical(datajags.region))';
% input.location = double(categorical(datajags.location))';
% input.study = double(categorical(datajags.sel))'; % sel?
% years = unique([double(categorical(datajags.sel)), datajags.midptyr], 'rows') ;
input.amtu5 = (datajags.amt_u5'); % str2double
input.amto5 = (datajags.amt_o5');

for i = 1:size(datajags,1)
if or(datajags.ages(i)=='1', or(datajags.ages(i)=='2', or(datajags.ages(i)=='3', datajags.ages(i)=='4')))
    input.offseta(i) = datajags.pop_time(i);
    input.offsetb(i) = datajags.pop_time(i);
    input.offsetc(i) = datajags.pop_time(i);
    input.offsetd(i) = datajags.pop_time(i);
elseif datajags.ages=='5'
    input.offseta(i) = datajags.popa(i)/(datajags.popa(i) + datajags.popb(i)).*datajags.pop_time(i);
    input.offsetb(i) = datajags.popb(i)/(datajags.popa(i) + datajags.popb(i)).*datajags.pop_time(i);
elseif datajags.ages=='6'
    input.offsetb(i) = datajags.popb(i)/(datajags.popb(i) + datajags.popc(i)).*datajags.pop_time(i);
    input.offsetc(i) = datajags.popc(i)/(datajags.popb(i) + datajags.popc(i)).*datajags.pop_time(i);
elseif datajags.ages=='7'
    input.offseta(i) = datajags.popa(i)/(datajags.popa(i) + datajags.popb(i) + datajags.popc(i)).*datajags.pop_time(i);
    input.offsetb(i) = datajags.popb(i)/(datajags.popa(i) + datajags.popb(i) + datajags.popc(i)).*datajags.pop_time(i);    
    input.offsetc(i) = datajags.popc(i)/(datajags.popa(i) + datajags.popb(i) + datajags.popc(i)).*datajags.pop_time(i);       
elseif datajags.ages=='8'
    input.offsetc(i) = datajags.popc(i)/(datajags.popc(i) + datajags.popd(i)).*datajags.pop_time(i);
    input.offsetd(i) = datajags.popd(i)/(datajags.popc(i) + datajags.popd(i)).*datajags.pop_time(i);    
else % if datajags.ages=='9'
    input.offseta(i) = datajags.popa(i)/100.*datajags.pop_time(i); 
    input.offsetb(i) = datajags.popb(i)/100.*datajags.pop_time(i); 
    input.offsetc(i) = datajags.popc(i)/100.*datajags.pop_time(i);
    input.offsetd(i) = datajags.popd(i)/100.*datajags.pop_time(i); 
end
end

input.I = size(datajags, 1); 
input.C = length(unique(data.region)); 
input.P = size(input.cov, 2);
input.W = eye(4);


