function EM = LyapunovExponent(x,d,p)
    N   = size(x,1);
    TOL = 0.1e-7;
    k   = 1:d;
    E   = zeros(N,1);
    for i = 1:N-d  % second quartile of data should be sufficiently evolved
        j         = 1:N-d;
        dist      = abs(x(i)-x(j));
        dist(i)   = inf;
        [~,indx]  = min(dist); % closest point!
         a         = ([(abs(x(i+k,:)-x(indx+k,:))>TOL).*(abs(x(i,:)-x(indx,:))>TOL)]==1);
        expn      = [log(abs(x(i+k,:)-x(indx+k,:)))-log(abs(x(i,:)-x(indx,:)))]./k';
        expn      = expn(a);
        if isempty(expn) == 1;expn = 0;end
        E(i) = sum(expn)/d;
    end
    EM = mean(E);
    if p == 1
        subplot(3,2,[1 2 3 4]);plot(x);yLabel('x');xlabel('Sample')
        subplot(3,2,[5 6]);plot(E);grid on
        xlim([1 N-d])
        ylim([0 max(E)+0.1*max(E)])
        xlabel('Exponent')
        ylabel('\lambda')
    end
end