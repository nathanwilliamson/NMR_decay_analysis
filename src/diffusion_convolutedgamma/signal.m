function [val, grad] = signal(b, mu, sigma, I0, Ib)

number_of_convolutions = numel(mu);

I = ones(size(b));

for current_convolution = 1:number_of_convolutions
    I = I .* ( mu(current_convolution) ./ (mu(current_convolution) + b * sigma(current_convolution)^2) ).^( mu(current_convolution)^2 / sigma(current_convolution)^2 );
end

val = I0 * I + Ib;

grad = nan(2 * number_of_convolutions + 2, numel(b));

for current_convolution = 1:number_of_convolutions
    grad(current_convolution, :) = I0 .* I .* mu(current_convolution) ./ sigma(current_convolution).^2 .* (2 * log(mu(current_convolution) ./ (mu(current_convolution) + b * sigma(current_convolution)^2)) + b * sigma(current_convolution).^2 ./ (mu(current_convolution) + b * sigma(current_convolution).^2));
    grad(number_of_convolutions + current_convolution, :) = - I0 * I .* mu(current_convolution).^2 ./ sigma(current_convolution).^2 .* (2 ./ sigma(current_convolution) .* log(mu(current_convolution) ./ (mu(current_convolution) + b * sigma(current_convolution)^2)) + 2 .* b .* sigma(current_convolution) ./ (mu(current_convolution) + b .* sigma(current_convolution)^2));
end

grad(2 * number_of_convolutions + 1, :) = I;
grad(2 * number_of_convolutions + 2, :) = 1;

end