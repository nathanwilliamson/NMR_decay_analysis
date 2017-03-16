function val = rand_mean_D(b)

% Dslope = (log(I(1))-log(I(2)))/(b(2)-b(1));
% val = 2 * Dslope * rand();

% Let b^-1 define the 'scale' of the problem when randomizing.
b = b(b > 0);
minD = 1 / b(end);
maxD = 1 / b(1);
val = minD + (maxD - minD) * rand();

end