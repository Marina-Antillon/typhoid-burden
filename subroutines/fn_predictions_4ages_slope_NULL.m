function [prelambda, lambda] = fn_predictions_4ages_slope_NULL(data, model, nsamples, covars, sample, burn, thin)
    
% data must be in wide format.

input.country = data.region; % double(categorical(data.region));
input.offseta = double(data.offseta);
input.offsetb = double(data.offsetb);
input.offsetc = double(data.offsetc);
input.offsetd = double(data.offsetd);

input.popa = double(data.popa);
input.popb = double(data.popb);
input.popc = double(data.popc);
input.popd = double(data.popd);

input.surv = double(data.surv==1);
input.amt_u5 = double(data.amt_u5);
input.amt_o5 = double(data.amt_o5);
input.part002 = double(data.part002);
input.part205 = double(data.part205);
input.part515 = double(data.part515);
input.part15o = double(data.part15o);

% predictions = zeros(size(data,1), 9, nsamples);
lambda = zeros(size(data,1), 9, nsamples);

    for i=1:nsamples
        
        if sample==1
            j=randsample(burn:thin:size(model.mu,1),1);
        else
            j=i;
        end

        par.alpha = model.alpha(j, :);
        par.zeta1 = model.zeta1(j, :); 
        par.zeta2 = model.zeta2(j, :); 
        par.zeta3 = model.zeta3(j, :); 
        
%          tmp = mvnrnd([model.mu(j,:); model.beta(j,:)'], [model.sigma1(j,:); model.sigma2(j,:); model.sigma3(j,:); model.sigma4(j,:)]); 
        
        par.mu = model.mu(j,:); % tmp(1); % 
        par.beta = model.beta(j,:); % tmp(2:4); % 

        par.betaprob = model.betaprob(j,:);
        par.betau5 = 1-exp(-par.betaprob.*input.amt_u5);
        par.betao5 = 1-exp(-par.betaprob.*input.amt_o5);
        par.psi = model.psi(j,:);

        surv_coef = ones(size(data,1), 1);
        surv_coef(input.surv==1) = par.psi;
        
        prelambda(input.country>0, 1, i) = exp(par.alpha(input.country(input.country>0))' + par.zeta1(input.country(input.country>0))');
        prelambda(input.country>0, 2, i) = exp(par.alpha(input.country(input.country>0))' + par.zeta2(input.country(input.country>0))');
        prelambda(input.country>0, 3, i) = exp(par.alpha(input.country(input.country>0))');
        prelambda(input.country>0, 4, i) = exp(par.alpha(input.country(input.country>0))' + par.zeta3(input.country(input.country>0))');
        
        lambda(input.country>0, 1, i) = exp(log(max(prelambda(input.country>0, 1, i), 1e-5)) + log(input.offseta(input.country>0)));
        lambda(input.country>0, 2, i) = exp(log(max(prelambda(input.country>0, 2, i), 1e-5)) + log(input.offsetb(input.country>0)));
        lambda(input.country>0, 3, i) = exp(log(max(prelambda(input.country>0, 3, i), 1e-5)) + log(input.offsetc(input.country>0)));
        lambda(input.country>0, 4, i) = exp(log(max(prelambda(input.country>0, 4, i), 1e-5)) + log(input.offsetd(input.country>0)));
         
        % Countries where there are no observations 
        prelambda(input.country==0, 1, i) = exp(par.mu + par.beta(1));
        prelambda(input.country==0, 2, i) = exp(par.mu + par.beta(2));
        prelambda(input.country==0, 3, i) = exp(par.mu);
        prelambda(input.country==0, 4, i) = exp(par.mu + par.beta(3));
        
        lambda(input.country==0, 1, i) = exp(log(max(prelambda(input.country==0, 1, i), 1e-5)) + log(input.offseta(input.country==0)));
        lambda(input.country==0, 2, i) = exp(log(max(prelambda(input.country==0, 2, i), 1e-5)) + log(input.offsetb(input.country==0)));
        lambda(input.country==0, 3, i) = exp(log(max(prelambda(input.country==0, 3, i), 1e-5)) + log(input.offsetc(input.country==0)));
        lambda(input.country==0, 4, i) = exp(log(max(prelambda(input.country==0, 4, i), 1e-5)) + log(input.offsetd(input.country==0)));
        
        if model.processes==2
        prelambda(:, 1, i) = prelambda(:, 1, i).*surv_coef.*par.betau5.*input.part002; % 'Young Children 0-1'
        prelambda(:, 2, i) = prelambda(:, 2, i).*surv_coef.*par.betau5.*input.part205; % 'Pre-K 2-4'
        prelambda(:, 3, i) = prelambda(:, 3, i).*surv_coef.*par.betao5.*input.part515; % 'School 5-14'
        prelambda(:, 4, i) = prelambda(:, 4, i).*surv_coef.*par.betao5.*input.part15o; % 'Adults 15+'
        end
        
        % Predict incidence after accounting for the observation process.
        if model.processes==2
        lambda(:, 1, i) = lambda(:, 1, i).*surv_coef.*par.betau5.*input.part002; % 'Young Children 0-1'
        lambda(:, 2, i) = lambda(:, 2, i).*surv_coef.*par.betau5.*input.part205; % 'Pre-K 2-4'
        lambda(:, 3, i) = lambda(:, 3, i).*surv_coef.*par.betao5.*input.part515; % 'School 5-14'
        lambda(:, 4, i) = lambda(:, 4, i).*surv_coef.*par.betao5.*input.part15o; % 'Adults 15+'
        end
        
        % rates of incidence per 100,000 people for combinations of age
        % groups. We weighted by the relative size of each of the age
        % subgroups within a combination.
        prelambda(:, 5, i) = prelambda(:, 1, i).*input.popa./(input.popa + input.popb) ...
            + prelambda(:, 2, i).*input.popb./(input.popa + input.popb); % 'Young & Pre-K 0-4'
        prelambda(:, 6, i) = prelambda(:, 2, i).*input.popb./(input.popb + input.popc) ...
            + prelambda(:, 3, i).*input.popc./(input.popb + input.popc) ; % 'Pre-K & School 2-14' 
        prelambda(:, 7, i) = prelambda(:, 1, i).*input.popa./(input.popa + input.popb + input.popc) ...
            + prelambda(:, 2, i).*input.popb./(input.popa + input.popb + input.popc) ...
            + prelambda(:, 3, i).*input.popc./(input.popa + input.popb + input.popc); % 'Young & Pre-K & School 0-14'
        prelambda(:, 8, i) = prelambda(:, 3, i).*input.popc./(input.popc + input.popd) ...
            + prelambda(:, 4, i).*input.popd./(input.popc + input.popd); % 'School & Adults 5+'       
        prelambda(:, 9, i) = prelambda(:, 1, i).*input.popa/100 + prelambda(:, 2, i).*input.popb/100 ...
            + prelambda(:, 3, i).*input.popc/100 + prelambda(:, 4, i).*input.popd/100; % 'all'
                  
        % because the sum of poisson processes is itself a Poisson process with the rate equal to the sum of the rates of the subprocesses.
        lambda(:, 5, i) = lambda(:, 1, i) + lambda(:, 2, i); % 'Young & Pre-K 0-4'
        lambda(:, 6, i) = lambda(:, 2, i) + lambda(:, 3, i); % 'Pre-K & School 2-14' 
        lambda(:, 7, i) = lambda(:, 1, i) + lambda(:, 2, i) + lambda(:, 3, i); % 'Young & Pre-K & School 0-14'
        lambda(:, 8, i) = lambda(:, 3, i) + lambda(:, 4, i); % 'School & Adults 5+'       
        lambda(:, 9, i) = lambda(:, 1, i) + lambda(:, 2, i) + lambda(:, 3, i) + lambda(:, 4, i); % 'all'

    end
    
end
