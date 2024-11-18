function plotTrajectory(Detect_Result)
    % plotTrajectory 绘制目标的极坐标和直角坐标动态轨迹，并保存为GIF
    % 输入：
    %   Detect_Result: 3xN的矩阵，第一行是目标距离，第三行是目标与y轴的角度
    
    % 获取数据
    distance = Detect_Result(1, :);         % 距离
    angle_xy = Detect_Result(3, :);         % 与y轴夹角
    angle_polar = -1 * angle_xy + 90;       % 转换为极坐标下的角度
    
    % 极坐标系动画
    filename_polar = 'detect_polar.gif';
    h = figure;
    for i = 1:length(distance)
        polarplot(angle_polar(i) * pi / 180, distance(i), 'bo');
        thetalim([0, 180]);
        % thetalim([75, 105]);
        % rlim([63500, 64500]);
        title('航迹动态显示');
        % thetaticklabels([]);
        % rticklabels([]);
        hold on;
        drawnow;
        
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind, cm] = rgb2ind(im, 256); 
        if i == 1 
            imwrite(imind, cm, filename_polar, 'gif', 'Loopcount', inf, 'DelayTime', 0.05); 
        else 
            imwrite(imind, cm, filename_polar, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05); 
        end
    end
    close(h);  % 关闭极坐标图窗口

    % 直角坐标系动画
    filename_xy = 'detect_xy.gif';
    X = distance .* sind(angle_xy);  % x 坐标
    Y = distance .* cosd(angle_xy);  % y 坐标
    b = figure;
    for i = 1:length(distance)
        plot(X(i), Y(i), 'bo');
        grid on;
        xlabel('水平距离（m）');
        ylabel('垂直距离（m）');
        axis([-4500 4500 0 70000]);  % 设置坐标轴范围
        title('航迹动态显示');
        hold on;
        drawnow;
        
        % Capture the plot as an image 
        frame = getframe(b); 
        im = frame2im(frame); 
        [imind, cm] = rgb2ind(im, 256); 
        if i == 1 
            imwrite(imind, cm, filename_xy, 'gif', 'Loopcount', inf, 'DelayTime', 0.05); 
        else 
            imwrite(imind, cm, filename_xy, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05); 
        end
    end
    close(b);  % 关闭直角坐标图窗口
end
