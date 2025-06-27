function beta0 = set_beta0_from_snr(SNR_dB, d_ref, alpha, sigma_w2)
% 根据目标SNR计算 beta0
% 输入:
%   SNR_dB     - 目标信噪比（单位：dB）
%   d_ref      - 参考距离，例如用户与UAV平均距离
%   alpha      - 路径损耗指数
%   sigma_w2   - 噪声功率
%
% 输出:
%   beta0      - 对应的参考信道增益

    SNR_linear = 10^(SNR_dB / 10);
    beta0 = SNR_linear * sigma_w2 * d_ref^alpha;
end