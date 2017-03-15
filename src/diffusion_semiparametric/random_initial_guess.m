function param_guess = random_initial_guess(model, lb, ub, b, I)

number_of_components = numel(model);

param_guess = nan(size(lb));
ind = 1;
for current_component = 1:number_of_components
    switch model{current_component}{1}
        case 'exponential'
            Drnd = randmeanD(b, I);
            Drnd = max(lb(ind), Drnd);
            Drnd = min(ub(ind), Drnd);
            param_guess(ind) = Drnd;
            ind = ind + 1;
        case 'stretchedexponential'
            Drnd = randmeanD(b, I);
            Drnd = max(lb(ind), Drnd);
            Drnd = min(ub(ind), Drnd);
            param_guess(ind) = Drnd;
            ind = ind + 1;
            
            betarnd = rand();
            betarnd = max(lb(ind), betarnd);
            betarnd = min(ub(ind), betarnd);
            param_guess(ind) = betarnd;
            ind = ind + 1;
        case 'lognormal'
            Mrnd = randmeanD(b,I);
            CVrnd = 0.05 + 0.95 * rand();
            
            murnd = log(Mrnd) - 1/2 * log(1 + CVrnd^2);
            murnd = max(lb(ind), murnd);
            murnd = min(ub(ind), murnd);
            param_guess(ind) = murnd;
            ind = ind + 1;
            
            sigmarnd = sqrt( log(1 + CVrnd^2) );
            sigmarnd = max(lb(ind), sigmarnd);
            sigmarnd = min(ub(ind), sigmarnd);
            param_guess(ind) = sigmarnd;
            ind = ind + 1;
        case 'gamma'
            Mrnd = randmeanD(b,I);
            CVrnd = 0.05 + 0.75 * rand();
            
            alpharnd = 1 / CVrnd^2;
            alpharnd = max(lb(ind), alpharnd);
            alpharnd = min(ub(ind), alpharnd);
            param_guess(ind) = alpharnd;
            ind = ind + 1;
            
            betarnd = 1 / (Mrnd * CVrnd^2);
            betarnd = max(lb(ind), betarnd);
            betarnd = min(ub(ind), betarnd);
            param_guess(ind) = betarnd;
            ind = ind + 1;
    end
end
baseline_rnd = abs(I(end)) * randn();
baseline_rnd = max(lb(ind), baseline_rnd);
baseline_rnd = min(ub(ind), baseline_rnd);
param_guess(ind) = baseline_rnd;
ind = ind + 1;

theta0 = rand(1, number_of_components);
theta0 = theta0 / sum(theta0); % This could be out of bounds, but will be corrected by the algorithm.
param_guess(ind:ind+number_of_components-1) = theta0;
param_guess(end) = 1.25 * max(I) * rand();

end

