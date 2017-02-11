% Blood culture sensitivity as a function of the volume of blood drawn from
% the patient.

clear
bc = readtable('../00-data/BC_sensitivity.csv');
% bc(end,:) = [];

repeats = {'Hoffman 1986', 'Gasem 1995'};
once = {'Gilman 1975', 'Guerra-Caceres 1979', 'Vallenas 1985', 'Gasem 2002', 'Wain 2008'}; 

input.obs2 = [27, 38-27, 23; 52, 1, 53];  
% I am assuming that for Wain 2008, in the end the sensivity for Blood 15
% mL was the same as for bone marrow (1-7/68)
input.all2 = sum(input.obs2, 2); 
input.vol2 = [3, 8; 3, 10]; 
input.study2 = [2; 5];

% Throw Vallenas in for the other analysis.
input.obs1 = [table2array(bc([1, 2, 3, 8, 9], {'positive'})); 59]; 
input.all1 = [sum(table2array(bc([1, 2, 3, 8, 9], {'positive', 'negative'})),2); 91]; 
input.vol1 = [table2array(bc([1, 2, 3, 8, 9], {'draw_midpt'})); 5]; 
input.study1 = [table2array(bc([1, 2, 3, 8, 9], {'study'})); 4];

input.I = length(input.obs1);
input.J = size(input.obs2,1);
input.K = length(unique(bc.study));

initial.alpha = 0.2;
initial.beta = 1;
initial.lambda= 0.2*ones(1, input.K);
% initial.lambda = 0.2;

nburnin = 5e3;
nsamples = 1e5;
nchains=1;
k=1; % thread number. In case I am running more than one jags process at a time.
doparallel = 0; % do not use parallelization
monitor = {'alpha', 'beta', 'lambda'}; % 'beta0', 'alpha',

fprintf( 'Running JAGS...\n' );
clock
tic
[~, ~, structArray_exp] = matjags( ...
    input, ...                     % Observed data   
    fullfile(pwd, 'sensitivity_estimate.jags'), ...    % File that contains model definition
    initial, k, ...                     % Initial values for latent variables
    'doparallel' , doparallel, ...      % Parallelization flag
    'nchains', nchains,...              % Number of MCMC chains
    'nburnin', nburnin,...              % Number of burnin steps
    'nsamples', nsamples, ...           % Number of samples to extract
    'thin', 10, ...                     % Thinning parameter
    'dic', 1, ...                       % Do the DIC?
    'monitorparams', monitor, ...       % List of latent variables to monitor
    'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
    'verbosity' , 2 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 0 );                    % clean up of temporary files?
toc


% Fit vs. Observed.

probability = bc.positive./(bc.positive+bc.negative);
errorlow = probability - betainv(0.025, bc.positive, bc.negative);
errorhigh = betainv(0.975, bc.positive, bc.negative) - probability;

ovol = bc.draw_midpt;
volval=unique(ovol);
for i=1:length(volval)
    tmp=find(ovol==volval(i));
    if length(tmp)>1
        for j=1:length(tmp)
            ovol(tmp(j)) = ovol(tmp(j)) + j*0.2;
        end
    end
end

% scatter(bc.volume, probability)
% errorbar(ovol, probability, errorlow, errorhigh, 'k.')
% Must calculate at various points, from 1-11

vol=[2:15];

clear prob
for i=1:10000
    rnd=randsample(1:length(structArray_exp.beta), 1);
%     lambda = gamrnd(structArray_exp.alpha(rnd), 1/structArray_exp.beta(rnd));
    lambda = lognrnd(structArray_exp.alpha(rnd), structArray_exp.beta(rnd)^(-0.5));
%     lambda = exp(structArray_exp.alpha(rnd));
    prob(i, :) = 1-exp(-lambda*vol);
    prob_mean = 1-exp(-exp(structArray_exp.alpha)*vol);
end

probsum = prctile(prob, [2.5 50 97.5]); 

upper = probsum(3, :);
lower = probsum(1, :);
med = probsum(2, :); 

probsum_mean = prctile(prob_mean, [2.5 50 97.5]); 

upper_mean = probsum_mean(3, :);
lower_mean = probsum_mean(1, :);
med_mean = probsum_mean(2, :); 

plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);

plot_variance(vol, lower, upper, [0.9 0.9 0.9])
hold on
plot_variance(vol, lower_mean, upper_mean, [0.7 0.7 0.7])
errorbar(ovol, probability, errorlow, errorhigh, 'k.')
axis([0 16 0 1])
ylabel('Blood Culture Sensitivity', 'Fontsize', 14)
xlabel('Sample Volume for Blood Culture (mL)', 'Fontsize', 14)
title({'Sample Volume vs Sensitivity'; 'All culture-confirmed cases'}, 'Fontsize', 14)
legend({'Population Response', 'Mean Prediction', 'Observed'}, 'Location', 'SouthOutside', 'Orientation', 'vertical', 'Fontsize', 14)  


months = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'}; 
today = clock;

% filename = strcat(['fit_obs_slope_re_' months{today(2)} sprintf('%2.0f', today(3))]);
    
% savefig(filename)

gamfit(lognrnd(structArray_exp.alpha, structArray_exp.beta.^(-0.5), 100000, 1))
gamfit(randsample(exp(structArray_exp.alpha), 100000, 1))

