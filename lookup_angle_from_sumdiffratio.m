function theta_closest = lookup_angle_from_sumdiffratio(AB_ybili, sum_diff_ratio, theta, theta_min, theta_max)
    % 检查 AB_ybili 和 theta 的长度是否一致
    if length(AB_ybili) ~= length(theta)
        error('AB_ybili 的长度必须与 theta 的长度一致');
    end

    % 限制查找范围为 [theta_min, theta_max]
    theta_range_idx = (theta >= theta_min) & (theta <= theta_max);
    theta_limited = theta(theta_range_idx);
    AB_ybili_limited = AB_ybili(theta_range_idx);

    % 计算 AB_ybili_limited 与 sum_diff_ratio 的差值
    diff = abs(AB_ybili_limited - sum_diff_ratio);

    % 找到差值最小的位置
    [~, idx] = min(diff);

    % 获取对应的 theta 值
    theta_closest = theta_limited(idx);
end
