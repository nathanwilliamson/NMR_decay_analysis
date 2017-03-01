function f = pdf_convolutedgamma(D,alpha,beta,nTerms)

beta1 = max(beta);
rho = sum(alpha);

logC = sum(alpha.*log(beta./beta1));

f = zeros(size(D));

gamma = nan(1,nTerms);

for k = 1:nTerms
    gamma_k = sum(alpha .* (1-beta/beta1).^k / k);
    gamma(k) = gamma_k;
end

delta = nan(1,nTerms+1);
delta0 = 1;
delta(0+1) = delta0;

for k = 1:nTerms
    delta_k = 1/k * sum( (1:k) .* gamma(1:k) .* delta(k:-1:1) );
    delta(k+1) = delta_k;
end

for k = 0:nTerms
    logterm = log(delta(k+1)) + (rho+k-1)*log(D) - beta1*D + (rho+k)*log(beta1) - gammaln(rho+k);
    term = exp(logterm);
    f = f + term;
end

f = f*exp(logC);

end
