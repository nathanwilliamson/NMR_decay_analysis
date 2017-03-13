function [val,grad] = signal_gamma(b,mu,sigma,I0,Ib)

I = (mu./(mu+b*sigma^2)).^(mu^2/sigma^2);

val = I0*I + Ib;

grad = nan(4,numel(b));

grad(1,:) = I0 .* I .* mu./sigma.^2.*(2*log(mu./(mu+b*sigma^2)) + b*sigma.^2./(mu+b*sigma.^2));
grad(2,:) = - I0 * I .* mu.^2./sigma.^2.*(2./sigma.*log(mu./(mu+b*sigma^2)) + 2.*b.*sigma./(mu+b.*sigma^2));
grad(3,:) = I;
grad(4,:) = ones(1,numel(b));

end