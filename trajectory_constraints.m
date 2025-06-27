function [c, ceq] = trajectory_constraints(qv, T, Vmax)
% qv: [2T×1] vector: [x1,...xT, y1,...yT]
x = qv(1:T);
y = qv(T+1:end);

c = zeros(T-1, 1);  % 非线性不等式约束
for t = 1:T-1
    dx = x(t+1) - x(t);
    dy = y(t+1) - y(t);
    c(t) = sqrt(dx^2 + dy^2) - Vmax;  % 限制每步最大速度
end

ceq = [];  % 无等式约束
end