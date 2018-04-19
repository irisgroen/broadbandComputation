% Hilbert of sum v sum of hilbert

dt = 0.001;
t = (dt:dt:1)';          % time
f = 1:10;                % frequency
a = ones(size(f));       % amplitude
p = rand(size(f))*2*pi;  % phase
%p = p*0;
a = rand(size(f));

sz = [length(t) length(f)];

T = bsxfun(@times, ones(sz), t);
F = bsxfun(@times, ones(sz), f);
P = bsxfun(@times, ones(sz), p);
A = bsxfun(@times, ones(sz), a);

X = A .* sin(2 * pi * T .* F + P);
x = mean(X, 2);

h = hilbert(mean(X,2));
H = hilbert(X);

v = mean(x.^2);
V = mean(X.^2);
figure(1),clf 

subplot(1,3,1)

plot(t, abs(h),'k',  t, -abs(h), 'k--', ...
    t, real(h), 'r', t, imag(h), 'g')

ylim(1.2*[-1 1])
title('Hilbert of mean')
legend('+abs', '-abs', 'real', 'imag')

subplot(1,3,2)

plot(t, mean(abs(H),2),'k',  t, -mean(abs(H),2), 'k--', ...
    t, mean(real(H),2), 'r', t, mean(imag(h),2), 'g')

ylim(1.2*[-1 1])
title('mean of Hilbert')
legend('+abs', '-abs', 'real', 'imag')

subplot(1,3,3)
plot(t, X, 'k'); hold on; plot(t, x, 'r', 'LineWidth', 3);
ylim(1.2*[-1 1])

