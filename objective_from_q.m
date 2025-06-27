function val = objective_from_q(qv, rho, chi, ak, bk, H, alpha, beta0, K, T, N, sigma_w2)
    x = qv(1:T); y = qv(T+1:end);
    eta = calc_eta(x, y, ak, bk, H, alpha, beta0, K, T);
    val = compute_total_rate(rho, chi, eta, sigma_w2, K, T, N);
end