n = 10;
repetitions = 100;
n_sig = 0;

for repetition = 1:repetitions
%     vector = complex(1, 1) + complex(randn(1, n), randn(1, n));
    vector = 1 +randn(1, n);
%     disp(sum(abs(vector)))
    p = signrank(vector);
    if p < 0.05
        n_sig = n_sig + 1;
    end
end

disp(n_sig)
