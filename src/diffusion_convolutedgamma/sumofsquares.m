function [val, grad] = sumofsquares(b, I, m, s, I0, Ib)

number_of_convolutions = numel(m);

[val_signal, grad_signal] = signal_convolutedgamma(b, m, s, I0, Ib);

val = sum( (I - val_signal).^2 );

grad = zeros(2 * number_of_convolutions + 2, 1);

for current_convolution = 1:number_of_convolutions
    grad(current_convolution, :) = - sum( 2 * (I - val_signal) .* grad_signal(current_convolution, :) );
    grad(number_of_convolutions + current_convolution, :) = - sum( 2 * (I - val_signal) .* grad_signal(number_of_convolutions + current_convolution, :) );
end

grad(2 * number_of_convolutions + 1, :) = -sum( 2 * (I - val_signal) .* grad_signal(2 * number_of_convolutions + 1, :) );
grad(2 * number_of_convolutions + 2, :) = -sum( 2 * (I - val_signal) .* grad_signal(2 * number_of_convolutions + 2, :) );

end

