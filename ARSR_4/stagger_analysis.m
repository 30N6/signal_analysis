pri_stagger_s = pri_stagger * 1e-6;

w_i = get_mti_weights(length(pri_stagger_s));
pri_cumsum = [0; cumsum(pri_stagger_s)];

f = linspace(0, 2e3, 1e4);
H_f = zeros(length(f), 1);

for i_f = 1:length(f)
    for ii = 1:length(pri_cumsum)
        H_f(i_f) = H_f(i_f) + w_i(ii) * exp(1j*2*pi*pri_cumsum(ii) * f(i_f));
    end
end

H = 10*log10(abs(H_f) / max(abs(H_f)));

figure(12);
plot(f, H);

function w_i = get_mti_weights(m)
    ii = 1:(m + 1);
    w_i = (-1).^(ii - 1) .* (factorial(m) ./ (factorial(m - ii + 1) .* factorial(ii - 1)));
end