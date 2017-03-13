function [val,grad] = sumofsquares_gamma(b,I,m,s,I0,Ib)

[val_signal,grad_signal] = signal_gamma(b,m,s,I0,Ib);

val = sum( (I - val_signal).^2 );

grad = zeros(4,1);

grad(1) = - sum( 2 * (I - val_signal) .* grad_signal(1,:) );
grad(2) = - sum( 2 * (I - val_signal) .* grad_signal(2,:) );
grad(3) = -sum( 2 * (I - val_signal) .* grad_signal(3,:) );
grad(4) = -sum( 2 * (I - val_signal) .* grad_signal(4,:) );

end

