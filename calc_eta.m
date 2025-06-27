function [eta,d_ave] = calc_eta(x, y, ak, bk, H, alpha, beta0, K, T)
    rng(27)
    eta = zeros(4, 4, K, T);
    sigma = 1e-0;
    K0 = 1000;
    for k = 1:K
        for t = 1:T
            d(k,t) = sqrt((x(t)-ak(k))^2 + (y(t)-bk(k))^2 + H^2);
            % disp(['d:',num2str(d)])
            for i = 1:4
                for j = 1:4
                    % h = sqrt(beta0 * (1/d)^alpha) * (sqrt(K0 / (K0 + 1)) + sqrt(1 / (K0 + 1)) / sqrt(2));
                    % h = sqrt(beta0 * (1/d)^alpha) * (sqrt(K0 / (K0 + 1)) + sqrt(1 / (K0 + 1)) * (randn + 1i*randn) / sqrt(2));abs(h)^2
                    % h = sqrt(beta0 * (1/d)^alpha) * sigma * (randn + 1i*randn)/sqrt(2);abs(h)^2
                    % h = sqrt(beta0 * (1/d)^alpha) * (1-1e-9 + 1e-9*rand)/sqrt(2);
                    h = sqrt(beta0 * (1/d(k,t))^alpha);%abs(h)^2
                    % disp(['eta:',num2str(abs(h)^2)])
                    eta(i,j,k,t) = abs(h)^2;
                end
            end
        end
    end
    d_ave = mean2(d);
end