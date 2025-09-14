function chi = optimize_all_chi(rho, chi, eta, sigma_w2, K, T, N)
    for k = 1:K
        for t = 1:T
            c0 = chi(:,k,t);
            obj = @(c) -compute_total_rate_single(rho(:,k,t), c, eta(:,:,k,t), sigma_w2);
            c = fmincon(obj, c0, [],[],[],[], zeros(N,1), ones(N,1), @(c) deal([],sum(c)-1), ...
                        optimoptions('fmincon', ...
                        'Algorithm','sqp', ...
                        'Display','off'));
            chi(:,k,t) = c;
        end
    end
end