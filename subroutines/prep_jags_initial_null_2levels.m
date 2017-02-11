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

formula = 'cases2 ~ 1 + ages + surv + (1 + ages|region)';   

options = statset('TolX', 1e-4);
try
init = fitglme(asd, formula, 'Distribution','Poisson','Link','log', 'Offset', log(asd.offset),'FitMethod','Laplace', 'DispersionFlag', false, 'CheckHessian', false, 'OptimizerOptions', options, 'EBOptions', options);
catch
init = fitglme(asd, formula, 'Distribution','Poisson','Link','log', 'Offset', log(asd.offset),'FitMethod','Laplace', 'DispersionFlag', false, 'CheckHessian', false, 'StartMethod', 'random', 'OptimizerOptions', options, 'EBOptions', options);
end

% parse so that I can give initial parameter values.
tmp = dataset2table(init.Coefficients);

% initial_null.mu = table2array(tmp(strcmp(tmp.Name, '(Intercept)'), 2));
initial_null.beta = repmat([table2array(tmp(strcmp(tmp.Name, '(Intercept)'), 2)), table2array(tmp(strcmp(tmp.Name, 'ages_1'), 2)), table2array(tmp(strcmp(tmp.Name, 'ages_2'), 2)), table2array(tmp(strcmp(tmp.Name, 'ages_4'), 2))], 3, 1);

initial_null.betaprob = 0.1747;
initial_null.psi = min(1, exp(table2array(tmp(strcmp(tmp.Name, 'surv'), 2)))); 

initial_null.gammapick = ones(1, size(input.cov, 2)); 

initial_null.gammaprior1 = zeros(1, size(input.cov, 2));
initial_null.gammaprior2 = zeros(1, size(input.cov, 2));
initial_null.gammaprior3 = zeros(1, size(input.cov, 2));
initial_null.gammaprior4 = zeros(1, size(input.cov, 2));

initial_null.taugamma = 1.1*ones(1, 4); 

%initial values related to the slope RE
initial_null.TauBraw = eye(4);
initial_null.Bcraw_r1(:,1) = initial_null.beta(1,1)*ones(input.C, 1); % this is the init for alpha
initial_null.Bcraw_r1(:,2) = initial_null.beta(1,2)*ones(input.C, 1);  % these are the init for zeta1, etc
initial_null.Bcraw_r1(:,3) = initial_null.beta(1,3)*ones(input.C, 1); 
initial_null.Bcraw_r1(:,4) = initial_null.beta(1,4)*ones(input.C, 1); 

initial_null.Bcraw_r2 = initial_null.Bcraw_r1;
initial_null.Bcraw_r3 = initial_null.Bcraw_r1;

% Get some plain average of point estimates of coefficients assign initial
% values to world hyperpriors.

initial_null.TauBeta = eye(4);
initial_null.omega = initial_null.beta(1,:);

