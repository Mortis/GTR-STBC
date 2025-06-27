function R = compute_total_rate(rho, chi, eta, sigma_w2, K, T, N)
    R = 0;
    for k = 1:K
        for t = 1:T
            num = 0;
            for i = 1:N
                for j = 1:N
                    % num = num + rho(i,k,t) * chi(j,k,t) * sqrt(eta(i,j,k,t));
                    num = num + rho(i,k,t) * chi(j,k,t) * eta(i,j,k,t);%
                    % num = num + rho(i,k,t) * chi(j,k,t) * sqrt(1e-6);
                end
            end
            num = num^2;

            % sum_rho2_chi2 = sum(rho(:,k,t).^2) * sum(chi(:,k,t).^2);
            sum_chi2 = sum(chi(:,k,t).^2);
            inner_eta = 0;
            for i2 = 1:N
                for j2 = 1:N
                    % inner_eta = inner_eta + rho(i2,k,t);
                    inner_eta = inner_eta + rho(i2,k,t) * eta(i2,j2,k,t);%
                end
            end
            denom = sigma_w2 * sum_chi2 * inner_eta;
            R = R + log2(1 + num / denom);
        end
    end
    % R = R / (K*T);
    R = R / T;
end