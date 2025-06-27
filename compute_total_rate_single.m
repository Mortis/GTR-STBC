function R = compute_total_rate_single(rho, chi, eta, sigma_w2)
    num = 0;
    for i = 1:4
        for j = 1:4
            num = num + rho(i) * chi(j) * eta(i,j);%
            % num = num + rho(i) * chi(j) * sqrt(eta(i,j));
            % num = num + rho(i) * chi(j) * sqrt(1e-6);
        end
    end
    num = num^2;

    % sum_rho2_chi2 = sum(rho.^2) * sum(chi.^2);
    sum_chi2 = sum(chi.^2);
    inner_eta = 0;
    for i2 = 1:4
        for j2 = 1:4
            inner_eta = inner_eta + rho(i2) * eta(i2,j2);%
            % inner_eta = inner_eta + rho(i2);
        end
    end
    denom = sigma_w2 * sum_chi2 * inner_eta;

    R = log2(1 + num / denom);
    
    % N = 4;
    % inner_eta = 0;
    % denom = 0;
    % for i2 = 1:N
    %     for j2 = 1:N
    %         inner_eta = inner_eta + rho(i2) * eta(i2,j2);
    %     end
    % end
    % num = 0;
    % for i = 1:N
    %     for j = 1:N
    %         num = num + rho(i) * chi(j) * sqrt(eta(i,j));
    %     end
    % end
    % num = num^2;
    % for i2 = 1:N
    %     for j2 = 1:N
    %         % denom = denom + rho(i2,k,t) * eta(i2,j2,k,t);
    %         denom = denom + rho(i).^2 * chi(j).^2 * sigma_w2 * sum(chi(:).^2) * inner_eta;
    %     end
    % end
    % % denom = sum_rho2_chi2 * sigma_w2 * sum_chi2 * inner_eta;
    % R = log2(1 + num / denom)
end