function [x_opt, y_opt] = optimize_trajectory(x0, y0, rho, chi, ak, bk, H, alpha, beta0, K, T, N, sigma_w2)
    % q0 = rand(2*T,1)+250;
    % x 线性从 0 到 100，y 固定
    % x0 = linspace(0, 100, T)';
    % y0 = ones(T, 1) * 50;
    Vmax = 25;
    q0 = [x0; y0];  % [x; y] 合并为 q0
    x0_init = 250;     % 起点x坐标
    y0_init = 250;     % 起点y坐标
    nonlcon = @(qv) trajectory_constraints_fixed(qv, T, Vmax, x0_init, y0_init);

    
options = optimoptions('fmincon', ...
   'Algorithm','sqp', ...
   'Display','off', ...
   'ConstraintTolerance', 1e-7, ...
   'StepTolerance', 1e-12, ...
   'OptimalityTolerance', 1e-7, ...
   'MaxIterations', 200);
    % options = optimoptions('fmincon','Algorithm', 'sqp','Display','off');
    % options = optimoptions('fmincon', ...
    % 'Algorithm', 'interior-point', ...
    % 'ConstraintTolerance', 1e-6, ...
    % 'StepTolerance', 1e-10, ...
    % 'MaxIterations', 100);

    obj = @(qv) -objective_from_q(qv, rho, chi, ak, bk, H, alpha, beta0, K, T, N, sigma_w2);
    % nonlcon = @(qv) trajectory_constraints(qv, T, Vmax);

    q_opt = fmincon(obj, q0, [], [], [], [], [], [], nonlcon, options);
    % q_opt = fmincon(obj, q0, [], [], [], [], [], [], [], options);
    x_opt = q_opt(1:T);
    y_opt = q_opt(T+1:end);
end