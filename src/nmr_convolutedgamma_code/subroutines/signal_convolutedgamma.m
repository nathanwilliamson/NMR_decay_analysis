function [val,grad] = signal_convolutedgamma(b,mu,sigma,I0,Ib)

nConv = numel(mu);

I = ones(size(b));

for i = 1:nConv
    I = I .* (mu(i)./(mu(i)+b*sigma(i)^2)).^(mu(i)^2/sigma(i)^2);
end

val = I0*I + Ib;

grad = nan(2*nConv+2,numel(b));

for i = 1:nConv
    grad(i,:) = I0 .* I .* mu(i)./sigma(i).^2.*(2*log(mu(i)./(mu(i)+b*sigma(i)^2)) + b*sigma(i).^2./(mu(i)+b*sigma(i).^2));
    grad(nConv+i,:) = - I0 * I .* mu(i).^2./sigma(i).^2.*(2./sigma(i).*log(mu(i)./(mu(i)+b*sigma(i)^2)) + 2.*b.*sigma(i)./(mu(i)+b.*sigma(i)^2));
end

grad(2*nConv+1,:) = I;
grad(2*nConv+2,:) = 1;

end