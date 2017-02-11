% Gather initial values for the parameters that will be fit via JAGS.
% These values are for the occasion when the random walk is initiated with
% a semi-saturated model (all predictors are included to model the 
% intercept, but none are included to model the slopes).

asd = data(or(data.ages=='1', or(data.ages=='2', or(data.ages=='3', data.ages=='4'))), :); 

asd.cases(asd.ages=='1') = round(asd.cases(asd.ages=='1')./asd.part002(asd.ages=='1'));
asd.cases(asd.ages=='2') = round(asd.cases(asd.ages=='2')./asd.part205(asd.ages=='2'));
asd.cases(asd.ages=='3') = round(asd.cases(asd.ages=='3')./asd.part515(asd.ages=='3'));
asd.cases(asd.ages=='4') = round(asd.cases(asd.ages=='4')./asd.part15o(asd.ages=='4'));

asd.cases2 = asd.cases*2;
asd = asd(~isnan(asd.cases),:);
asd.surv = double(asd.survtype1=='P'); % table2array(datajags(:,{'surveillance1','surveillance2'}));

asd.ages = mergecats(asd.ages, {'4', '5'}, '4');
asd.ages = mergecats(asd.ages, {'4', '6'}, '4');
asd.ages = mergecats(asd.ages, {'4', '7'}, '4');
asd.ages = mergecats(asd.ages, {'4', '8'}, '4');
asd.ages = mergecats(asd.ages, {'4', '9'}, '4');

% *Tell Matlab to make the 5-15 group the referrent.*
asd.ages = reordercats(asd.ages, {'3', '1', '2', '4'}); 
     
chosen = 1:size(input.cov, 2);
covar_mod = covar{chosen(1)}; 
    for i = 2:length(chosen)
        covar_mod = strcat(covar_mod, '+', covar{chosen(i)});
    end
   
formula = strcat('cases2 ~ 1 + ages + surv + ', covar_mod, ' + ', ' (1 + ages|region)'); % , covar_int,

options = statset('TolX', 1e-4);
try
init = fitglme(asd, formula, 'Distribution','Poisson','Link','log', 'Offset', log(asd.offset),'FitMethod','Laplace', 'DispersionFlag', false, 'CheckHessian', false, 'OptimizerOptions', options, 'EBOptions', options);
catch
init = fitglme(asd, formula, 'Distribution','Poisson','Link','log', 'Offset', log(asd.offset),'FitMethod','Laplace', 'DispersionFlag', false, 'CheckHessian', false, 'StartMethod', 'random', 'OptimizerOptions', options, 'EBOptions', options);
end

% parse so that I can give gamma_int1 and gamma_int2 initial values.
tmp = dataset2table(init.Coefficients);
    cov_keep = (xor(xor(xor(strcmp(tmp.Name, {'(Intercept)'}), strcmp(tmp.Name,{'ages_1'})), strcmp(tmp.Name,{'ages_2'})), strcmp(tmp.Name,{'ages_4'}))~=1);
    cov_names = tmp.Name(cov_keep); 


initial_sat.mu = table2array(tmp(strcmp(tmp.Name, '(Intercept)'), 2));
initial_sat.beta = [table2array(tmp(strcmp(tmp.Name, 'ages_1'), 2)), ...
                table2array(tmp(strcmp(tmp.Name, 'ages_2'), 2)), ...
                table2array(tmp(strcmp(tmp.Name, 'ages_4'), 2))] ;
initial_sat.betaprob = 0.1747;
initial_sat.psi = min(1, exp(table2array(tmp(strcmp(tmp.Name, 'surv'), 2)))); % exp(priors(1:2))';

initial_sat.gammapick = 2*ones(1, size(input.cov, 2)); 

initial_sat.gammaprior1 = zeros(1, size(input.cov, 2));
initial_sat.gammaprior2 = zeros(1, size(input.cov, 2));
initial_sat.gammaprior3 = zeros(1, size(input.cov, 2));
initial_sat.gammaprior4 = zeros(1, size(input.cov, 2));

for i=1:length(chosen) 
   initial_sat.gammaprior1(1, i) = table2array(tmp(strcmp(tmp.Name, covar{chosen(i)}), 2)) ;
end

initial_sat.gammaprior1 = max(min(initial_sat.gammaprior1, 3.9), -3.9); 
initial_sat.taugamma = [3.99, 2, 2, 2]; 

%initial values related to the slope RE
initial_sat.TauBraw = eye(4);
initial_sat.Bcraw(:,1) = initial_sat.mu*ones(input.C, 1); 
initial_sat.Bcraw(:,2) = initial_sat.beta(1)*ones(input.C, 1);  
initial_sat.Bcraw(:,3) = initial_sat.beta(2)*ones(input.C, 1); 
initial_sat.Bcraw(:,4) = initial_sat.beta(3)*ones(input.C, 1); 
