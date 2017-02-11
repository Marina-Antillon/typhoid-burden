function [prelambda, lambda] = fn_predictions_4ages_obs(data, model, nsamples, covars, sample, burn, thin)

input.surv = double(data.surv==1);

input.offseta = double(data.offseta);
input.offsetb = double(data.offsetb);
input.offsetc = double(data.offsetc);
input.offsetd = double(data.offsetd);
input.offsete = double(data.offsete);
input.offsetf = double(data.offsetf);
input.offsetg = double(data.offsetg);
input.offseth = double(data.offseth);
input.offseti = double(data.offseti);

input.amt_u5 = double(data.amt_u5);
input.amt_o5 = double(data.amt_o5);
input.part002 = double(data.part002);
input.part205 = double(data.part205);
input.part515 = double(data.part515);
input.part15o = double(data.part15o);

prelambda = zeros(size(data,1), 9, nsamples);
lambda = zeros(size(data,1), 9, nsamples);

    for i=1:nsamples
        
        if sample==1
            j=randsample(burn:thin:size(model.mu,1),1);
        else
            j=i;
        end
        
        % Get a draw of the incidence parameter from a gamma distribution
        % parameterized by the cases and the population 
        % + 1 so that observations=0 have an adequate CI (by the Rule of
        % Three)

        prelambda(:, 1, i) = gamrnd(data.casesa+1, 1./(data.offseta+1));
        prelambda(:, 2, i) = gamrnd(data.casesb+1, 1./(data.offsetb+1));
        prelambda(:, 3, i) = gamrnd(data.casesc+1, 1./(data.offsetc+1));
        prelambda(:, 4, i) = gamrnd(data.casesd+1, 1./(data.offsetd+1));
        prelambda(:, 5, i) = gamrnd(data.casese+1, 1./(data.offsete+1));
        prelambda(:, 6, i) = gamrnd(data.casesf+1, 1./(data.offsetf+1));
        prelambda(:, 7, i) = gamrnd(data.casesg+1, 1./(data.offsetg+1));
        prelambda(:, 8, i) = gamrnd(data.casesh+1, 1./(data.offseth+1));
        prelambda(:, 9, i) = gamrnd(data.casesi+1, 1./(data.offseti+1));
        
        lambda(:, 1, i) = prelambda(:, 1, i).*input.offseta; 
        lambda(:, 2, i) = prelambda(:, 2, i).*input.offsetb; 
        lambda(:, 3, i) = prelambda(:, 3, i).*input.offsetc; 
        lambda(:, 4, i) = prelambda(:, 4, i).*input.offsetd; 
        lambda(:, 5, i) = prelambda(:, 5, i).*input.offsete; 
        lambda(:, 6, i) = prelambda(:, 6, i).*input.offsetf; 
        lambda(:, 7, i) = prelambda(:, 7, i).*input.offsetg; 
        lambda(:, 8, i) = prelambda(:, 8, i).*input.offseth; 
        lambda(:, 9, i) = prelambda(:, 9, i).*input.offseti; 

        % Adjust for passive/active surveillance (using the adjusted model's estimated surveillance effect)
        par.psi = model.psi(j,:);
        % Adjust for blood culture sensitivity
        par.betaprob = model.betaprob(j,:);
        par.betau5 = 1-exp(-par.betaprob.*input.amt_u5);
        par.betao5 = 1-exp(-par.betaprob.*input.amt_o5);        
        
        surv_coef = ones(size(data,1), 1);
        surv_coef(input.surv==1) = par.psi;
        
        % Now adjust for participation rate, and the parameters for blood
        % culture sensitivity and surveillance coefficient.
        if model.processes==1; 
        prelambda(:, 1, i) = prelambda(:, 1, i)./surv_coef./par.betau5./input.part002; % 'Young Children 0-1'
        prelambda(:, 2, i) = prelambda(:, 2, i)./surv_coef./par.betau5./input.part205; % 'Pre-K 2-4'
        prelambda(:, 3, i) = prelambda(:, 3, i)./surv_coef./par.betao5./input.part515; % 'School 5-14'
        prelambda(:, 4, i) = prelambda(:, 4, i)./surv_coef./par.betao5./input.part15o; % 'Adults 15+'

        prelambda(:, 5, i) = prelambda(:, 5, i)./surv_coef./par.betao5./input.part515; % 'Young & Pre-K 0-4'
        prelambda(:, 6, i) = prelambda(:, 6, i)./surv_coef./par.betao5./input.part515; % 'Pre-K & School 2-14'
        prelambda(:, 7, i) = prelambda(:, 7, i)./surv_coef./par.betao5./input.part515; % 'Young & Pre-K & School 0-14'
        prelambda(:, 8, i) = prelambda(:, 8, i)./surv_coef./par.betao5./input.part515; % 'School & Adults 5+'
        prelambda(:, 9, i) = prelambda(:, 9, i)./surv_coef./par.betao5./input.part515; % 'all'
        
        lambda(:, 1, i) = lambda(:, 1, i)./surv_coef./par.betau5./input.part002; % 'Young Children 0-1'
        lambda(:, 2, i) = lambda(:, 2, i)./surv_coef./par.betau5./input.part205; % 'Pre-K 2-4'
        lambda(:, 3, i) = lambda(:, 3, i)./surv_coef./par.betao5./input.part515; % 'School 5-14'
        lambda(:, 4, i) = lambda(:, 4, i)./surv_coef./par.betao5./input.part15o; % 'Adults 15+'

        lambda(:, 5, i) = lambda(:, 5, i)./surv_coef./par.betao5./input.part515; % 'Young & Pre-K 0-4'
        lambda(:, 6, i) = lambda(:, 6, i)./surv_coef./par.betao5./input.part515; % 'Pre-K & School 2-14'
        lambda(:, 7, i) = lambda(:, 7, i)./surv_coef./par.betao5./input.part515; % 'Young & Pre-K & School 0-14'
        lambda(:, 8, i) = lambda(:, 8, i)./surv_coef./par.betao5./input.part515; % 'School & Adults 5+'
        lambda(:, 9, i) = lambda(:, 9, i)./surv_coef./par.betao5./input.part515; % 'all'
        end
        
    end
    
end