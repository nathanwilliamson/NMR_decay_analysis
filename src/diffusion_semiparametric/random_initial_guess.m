function param_guess = random_initial_guess( input_args )

param_guess = [];
ind = 1;
for currentComponent = 1:number_of_components
    switch model{currentComponent}{1}
        case 'exponential'
            Drnd = randmeanD(b,I);
            Drnd = max(lb(ind),Drnd);
            Drnd = min(ub(ind),Drnd);
            param_guess(ind) = Drnd;
            ind = ind + 1;
        case 'stretchedexponential'
            Drnd = randmeanD(b,I);
            Drnd = max(lb(ind),Drnd);
            Drnd = min(ub(ind),Drnd);
            param_guess(ind) = Drnd;
            ind = ind + 1;
            
            betarnd = rand();
            betarnd = max(lb(ind),betarnd);
            betarnd = min(ub(ind),betarnd);
            param_guess(ind) = betarnd;
            ind = ind + 1;
        case 'lognormal'
            Mrnd = randmeanD(b,I);
            CVrnd = 0.05 + 1 * rand();
            
            murnd = log(Mrnd) - 1/2*log(1+CVrnd^2);
            murnd = max(lb(ind),murnd);
            murnd = min(ub(ind),murnd);
            param_guess(ind) = murnd;
            ind = ind + 1;
            
            sigmarnd = sqrt(log(1+CVrnd^2));
            sigmarnd = max(lb(ind),sigmarnd);
            sigmarnd = min(ub(ind),sigmarnd);
            param_guess(ind) = sigmarnd;
            ind = ind + 1;
        case 'gamma'
            Mrnd = randmeanD(b,I);
            CVrnd = 0.05 + 0.75 * rand();
            
            alpharnd = 1/CVrnd^2;
            alpharnd = max(lb(ind),alpharnd);
            alpharnd = min(ub(ind),alpharnd);
            param_guess(ind) = alpharnd;
            ind = ind + 1;
            
            betarnd = 1/(Mrnd*CVrnd^2);
            betarnd = max(lb(ind),betarnd);
            betarnd = min(ub(ind),betarnd);
            param_guess(ind) = betarnd;
            ind = ind + 1;
    end
end
brnd = abs(I(end))*randn();
brnd = max(lb(ind),brnd);
brnd = min(ub(ind),brnd);
param_guess(ind) = brnd;

theta0 = rand(1,number_of_components);
theta0 = theta0 / sum(theta0); % This can be out of bounds, but will be corrected by the algorithm.
param_guess = [param_guess theta0];
param_guess = [param_guess max(I)*rand()];


end

