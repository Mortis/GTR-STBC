function rho = optimize_all_rho(rho, chi, eta, sigma_w2, K, T, N)
    for k = 1:K
        for t = 1:T
            r0 = rho(:,k,t);
            obj = @(r) -compute_total_rate_single(r, chi(:,k,t), eta(:,:,k,t), sigma_w2);
            r = fmincon(obj, r0, [],[],[],[], zeros(N,1), ones(N,1), @(r) deal([],sum(r)-1), ...
                        optimoptions('fmincon','Algorithm','sqp','Display','off'));
            rho(:,k,t) = r;
        end
    end
end

% function [rho, all_fval_hist] = optimize_all_rho(rho, chi, eta, sigma_w2, K, T, N)
%     all_fval_hist = cell(K,T);  % 存储每个用户、每个时隙的 fval 轨迹
% 
%     for k = 1:K
%         for t = 1:T
%             r0 = rho(:,k,t);
%             obj = @(r) -compute_total_rate_single(r, chi(:,k,t), eta(:,:,k,t), sigma_w2);
% 
%             % 清空上一次的历史记录（由 OutputFcn 写入 base 工作区）
%             evalin('base', 'clear last_fval_history');
% 
%             % 设置 fmincon options，绑定输出函数
%             options = optimoptions('fmincon', ...
%                 'Display','iter', ...
%                 'OutputFcn', @store_obj_history);
% 
%             % 执行优化
%             r = fmincon(obj, r0, [],[],[],[], zeros(N,1), ones(N,1), @(r) deal([],sum(r)-1), options);
%             rho(:,k,t) = r;
% 
%             % 获取该次优化的目标函数值轨迹
%             fval_hist = evalin('base', 'last_fval_history');
%             all_fval_hist{k,t} = fval_hist;
%         end
%     end
% end
