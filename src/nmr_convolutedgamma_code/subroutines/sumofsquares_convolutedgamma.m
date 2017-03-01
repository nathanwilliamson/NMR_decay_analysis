function [val,grad] = sumofsquares_convolutedgamma(b,I,m,s,I0,Ib)

nConv = numel(m);

[val_signal,grad_signal] = signal_convolutedgamma(b,m,s,I0,Ib);

val = sum( (I - val_signal).^2 );

grad = zeros(2*nConv+2,1);

for i = 1:nConv
    grad(i,:) = - sum( 2 * (I - val_signal) .* grad_signal(i,:) );
    grad(nConv+i,:) = - sum( 2 * (I - val_signal) .* grad_signal(nConv+i,:) );
end

grad(2*nConv+1,:) = -sum( 2 * (I - val_signal) .* grad_signal(2*nConv+1,:) );
grad(2*nConv+2,:) = -sum( 2 * (I - val_signal) .* grad_signal(2*nConv+2,:) );

end

