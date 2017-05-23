function y = signal_bandwidth(x,fs,sc)
    % Bandwidth function

    N1 = 1024;    % Points in FFT
    N2 = N1/2;    % Half FFT length
    fx = fft(x,N1);  % FFT of x
    f = abs(fx(1:N2)); % Find magnitude spectrum
    magSq = f.^2;     % Magnitude squared

    for i=1:N2
        ksc = ((i/N1*fs) - sc)^2;   % (k-sc) squared
        num(i) = ksc * magSq(i);    % Multiply by magnitude squared      
    end
    num = sum(num);
    den = sum(magSq);

    y = sqrt(num/den);
end