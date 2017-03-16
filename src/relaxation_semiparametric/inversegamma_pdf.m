function val = inversegamma_pdf(x, alpha, beta)

val = beta^alpha / gamma(alpha) * x.^(-alpha-1) .* exp(-beta ./ x);
val(x == 0) = 0;

end

