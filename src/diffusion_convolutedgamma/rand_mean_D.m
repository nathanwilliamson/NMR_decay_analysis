function val = rand_mean_D(b)

b = b(b > 0);
minD = 1 / b(end);
maxD = 1 / b(1);
val = minD + (maxD - minD) * rand();

end