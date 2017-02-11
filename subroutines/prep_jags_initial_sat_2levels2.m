% Gather initial values for the parameters that will be fit via JAGS.
% These values are for the occasion when the random walk is initiated with
% a NULL model.

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
% covar_int = strcat(strcat(' ages:', covar{chosen(1)}), ' + '); 
% for i = 2:length(chosen)
%     covar_int = strcat(covar_int, strcat(' ages:', covar{chosen(i)}), ' + ');  % what if there is only one covariate?  
% end    
formula = strcat('cases2 ~ 1 + ages + surv + ', covar_mod, ' + ', ' (1 + ages|region)'); % , covar_int,

options = statset('TolX', 1e-4);
try
init = fitglme(asd, formula, 'Distribution','Poisson','Link','log', 'Offset', log(asd.offset),'FitMethod','Laplace', 'DispersionFlag', false, 'CheckHessian', false, 'OptimizerOptions', options, 'EBOptions', options);
catch
init = fitglme(asd, formula, 'Distribution','Poisson','Link','log', 'Offset', log(asd.offset),'FitMethod','Laplace', 'DispersionFlag', false, 'CheckHessian', false, 'StartMethod', 'random', 'OptimizerOptions', options, 'EBOptions', options);
end

% parse so that I can give initial parameter values.
tmp = dataset2table(init.Coefficients);

    cov_keep = (xor(xor(xor(strcmp(tmp.Name, {'(Intercept)'}), strcmp(tmp.Name,{'ages_1'})), strcmp(tmp.Name,{'ages_2'})), strcmp(tmp.Name,{'ages_4'}))~=1);
%     cov_keep((end-length(covar)):end,:) = 0;
    cov_names = tmp.Name(cov_keep); 

% initial_sat.mu = table2array(tmp(strcmp(tmp.Name, '(Intercept)'), 2));
initial_sat.beta = repmat([table2array(tmp(strcmp(tmp.Name, '(Intercept)'), 2)), table2array(tmp(strcmp(tmp.Name, 'ages_1'), 2)), table2array(tmp(strcmp(tmp.Name, 'ages_2'), 2)), table2array(tmp(strcmp(tmp.Name, 'ages_4'), 2))], 6, 1);

initial_sat.betaprob = 0.1747;
initial_sat.psi = min(1, exp(table2array(tmp(strcmp(tmp.Name, 'surv'), 2)))); 

initial_sat.gammapick = ones(1, size(input.cov, 2)); 

initial_sat.gammaprior1 = zeros(1, size(input.cov, 2));
initial_sat.gammaprior2 = zeros(1, size(input.cov, 2));
initial_sat.gammaprior3 = zeros(1, size(input.cov, 2));
initial_sat.gammaprior4 = zeros(1, size(input.cov, 2));

for i=1:length(chosen) 
   initial_sat.gammaprior1(1, i) = table2array(tmp(strcmp(tmp.Name, covar{chosen(i)}), 2)) ;
%    initial.gammaprior2(1, i) = table2array(tmp_int1(strcmp(tmp_int1.Name, covar{chosen(i)}), 2)) ;
%    initial.gammaprior3(1, i) = table2array(tmp_int2(strcmp(tmp_int2.Name, covar{chosen(i)}), 2)) ;
%    initial.gammaprior4(1, i) = table2array(tmp_int3(strcmp(tmp_int3.Name, covar{chosen(i)}), 2)) ;
end

initial_sat.gammaprior1 = max(min(initial_sat.gammaprior1, 3.9), -3.9); 
% initial.gammaprior2 = max(min(initial.gammaprior2, 3.9), -3.9); 
% initial.gammaprior3 = max(min(initial.gammaprior3, 3.9), -3.9); 
% initial.gammaprior4 = max(min(initial.gammaprior4, 3.9), -3.9); 

% initial.mugamma = zeros(1,4); 
initial_sat.taugamma = [3.99, 2, 2, 2];

%initial values related to the slope RE
initial_sat.TauBraw = eye(4);
initial_sat.Bcraw_r1(:,1) = initial_sat.beta(1,1)*ones(input.C, 1); % this is the init for alpha
initial_sat.Bcraw_r1(:,2) = initial_sat.beta(1,2)*ones(input.C, 1);  % these are the init for zeta1, etc
initial_sat.Bcraw_r1(:,3) = initial_sat.beta(1,3)*ones(input.C, 1); 
initial_sat.Bcraw_r1(:,4) = initial_sat.beta(1,4)*ones(input.C, 1); 

initial_sat.Bcraw_r2 = initial_sat.Bcraw_r1;
initial_sat.Bcraw_r3 = initial_sat.Bcraw_r1;
initial_sat.Bcraw_r4 = initial_sat.Bcraw_r1;
initial_sat.Bcraw_r5 = initial_sat.Bcraw_r1;
initial_sat.Bcraw_r6 = initial_sat.Bcraw_r1;

% Get some plain average of point estimates of coefficients assign initial
% values to world hyperpriors.

initial_sat.TauBeta = eye(4);
initial_sat.omega = initial_sat.beta(1,:);

