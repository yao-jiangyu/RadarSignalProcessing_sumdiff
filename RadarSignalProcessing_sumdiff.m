% 雷达测速测距测角
clc;clear;close all;

%% 参数设置
% 基础参数
c = 3.0e8; % 光速(m/s)
Fc = 35e9; % 雷达射频
Br = 10e6; % 发射信号带宽
fs = 20*1e6; % 采样频率
PRF = 2e3; % 脉冲重复频率
PRT = 1/PRF; % 脉冲重复周期
lamda = c/Fc; % 雷达工作波长，用于计算多普勒频移
N_pulse = 128; % 回波脉冲数
N_sample = round(fs*PRT); % 每个脉冲周期的采样点数；
Tr = 3*1e-6; % 发射信号时宽
t1 = (0:1/fs:(N_sample-1)/fs); % 时间序列
RangeMax = c*t1(end)/2;% 最大不模糊距离
Range = c*t1/2; % 距离序列
Vmax = lamda*PRF/2; % 最大可检测速度
Velovity = -Vmax/2:Vmax/N_pulse:Vmax/2-Vmax/N_pulse; % 速度序列
searching_doa = -15:0.01:15; % 角度搜索区间

% 阵列参数
M = 16; % 阵元数量
SourceNum = 1; % 信号源数量
d = lamda/2; % 阵元间隔
d_LinearArray = (0:M-1)'*d; % 阵元间距
SNR = 10;
SNR = 10^(SNR/10); % 信杂比 Pt=10w RCS=20m^2

% 场景参数
V = 50;
T = 1; % 采样间隔，没有采用PRT是因为：T=PRT时，nT=300001,跑一次代码太久
nT = size((-4e3:V*T:4e3),2); % 采样帧数
Xk = [-4e3:V*T:4e3; 64e3*ones(1,nT); V*ones(1,nT); zeros(1,nT)]; % 目标状态变化
for i = 1:nT
    r = sqrt(Xk(1,i)^2 + Xk(2,i)^2); % 径向距离
    v = -(Xk(1,i) * Xk(3,i) + Xk(2,i) * Xk(4,i)) / r; % 径向速度
    phi = -(atan2d(Xk(2,i), Xk(1,i))-90); % 角度
    Zk(:,i) = [r;v;phi]; % i时刻目标极坐标状态
end

%% 波束鉴角曲线的生成
theta = -90:0.01:90;
theta1 = -3;             %波形A指向的方向（度）
theta2 = 3;             %波束B指向的方向
theta_min = -3.8;
theta_max = 3.8;
look_a = exp(1j*2*pi*d_LinearArray*sind(theta)/lamda);    %导向矢量
w_1 = exp(1j*2*pi*d_LinearArray*sind(theta1)/lamda); %波束A加权权向量
w_2 = exp(1j*2*pi*d_LinearArray*sind(theta2)/lamda); %波束B加权权向量
yA = abs(w_1'*look_a);                                %波束A的方向图
yB = abs(w_2'*look_a);                                %波束B的方向图
ABSum = yA+yB;                                    %和波束的方向图
ABDiff = yA-yB;                                   %差波束的方向图
AB_ybili = ABDiff./ABSum;                              % 差和比
% 绘制两波束
figure;
plot(theta,(yA/max(yA)),'linewidth',1);   %绘制波束A
hold on;
plot(theta,(yB/max(yB)),'linewidth',1);   %绘制波束B
xlabel('方位角/°');
ylabel('归一化方向图');
legend('波束A','波束B');
title('波束A、B示意图');
axis tight;
grid on;
% 绘制和差波束
figure;
plot(theta,ABSum,'linewidth',1);   %绘制和波束
hold on;
plot(theta,ABDiff,'linewidth',1);   %绘制差波束
xlabel('方位角/°');
ylabel('功率增益');
legend('和波束','差波束');
title('和差波束示意图');
axis tight;
grid on;
% 绘制鉴角曲线
figure
plot(theta,AB_ybili);
xlim([theta_min theta_max]);
xlabel('方位角/°');
ylabel('差和比');
title('鉴角曲线');
grid on;

%% 主程序

%初始化数组
Detect_Result = zeros(3,nT); % 最终测量结果

signal_LFM = zeros(M,N_pulse,N_sample); % 信号矩阵

signal_i = ones(N_pulse,N_sample); % 中间累加矩阵变量
y1_out = ones(N_pulse,N_sample);
y2_out = ones(N_pulse,N_sample);

FFT_y1out_all = ones(N_pulse,N_sample,nT); % 保存每个nT时刻下MTD、CFAR结果
FFT_y2out_all = ones(N_pulse,N_sample,nT);
RDM_mask_A_all = ones(N_pulse,N_sample,nT);  
RDM_mask_B_all = ones(N_pulse,N_sample,nT);

%匹配滤波系数生成
sr = rectpuls(t1-Tr/2,Tr).*exp(1j*pi*(Br/Tr).*(t1-Tr/2).^2);%LFM发射信号
win = hamming(N_sample)'; %匹配滤波加窗
win2 = repmat(hamming(N_pulse),1,N_sample);  %MTD加窗
h_w = fliplr(conj(sr)).*win;
h_w_freq = fft(h_w);

%噪声
clutter = sqrt(2)/2*randn(M,N_sample)+sqrt(2)/2*1i*randn(M,N_sample); % 噪声

for t = 1:nT
    data = Zk(:,t); % 读取目标真实位置
    a_tar_LinearArray = exp(1j*2*pi*d_LinearArray*sind(data(3))/lamda); % 期望信号的导向矢量，线性阵列

    for i_n = 1:N_pulse
        ta = (i_n-1)*PRT;
        tao = 2*(data(1)-data(2).*(ta+t1))/c;
        signal_i(i_n,:) = SNR.*rectpuls(t1-tao-Tr/2,Tr).*exp(1j*2*pi*Fc*(t1-tao-Tr/2)+1j*pi*(Br/Tr).*(t1-tao-Tr/2).^2);

        signal_LFM(:,i_n,:) = a_tar_LinearArray * signal_i(i_n,:) + clutter;
        st = squeeze(signal_LFM(:,i_n,:));

        y1 = w_1'*st;                                %波束A回波
        y2 = w_2'*st;                                %波束B回波

        % 脉冲压缩
        y1_out(i_n,:) = ifft(fft(y1,N_sample,2).*h_w_freq,N_sample,2);
        y2_out(i_n,:) = ifft(fft(y2,N_sample,2).*h_w_freq,N_sample,2);

    end

    %% MTD
    FFT_y1out = fftshift(fft(y1_out.*win2),1);
    FFT_y2out = fftshift(fft(y2_out.*win2),1);
    FFT_y1out_all(:,:,t) = FFT_y1out;
    FFT_y2out_all(:,:,t) = FFT_y2out;

    % figure
    % mesh(abs(FFT_y1out))
    % mesh(Range,Velovity,abs(FFT_y1out))
    % figure
    % mesh(abs(FFT_y2out))
    % mesh(Range,Velovity,abs(FFT_y2out))
    %% CA-CFAR

    numGuard = 2; % # of guard cells
    numTrain = numGuard*2; % # of training cells
    P_fa = 1e-5; % desired false alarm rate 
    SNR_OFFSET = -5; % dB

    RDM_dB_y1 = 10*log10(abs(FFT_y1out)/max(max(abs(FFT_y1out))));
    RDM_dB_y2 = 10*log10(abs(FFT_y2out)/max(max(abs(FFT_y2out))));

    % 对波束 A 和波束 B 分别执行 CA-CFAR 检测
    [RDM_mask_A, cfar_ranges_A, cfar_dopps_A, K_A] = ca_cfar(RDM_dB_y1, numGuard, numTrain, P_fa, SNR_OFFSET);
    [RDM_mask_B, cfar_ranges_B, cfar_dopps_B, K_B] = ca_cfar(RDM_dB_y2, numGuard, numTrain, P_fa, SNR_OFFSET);
    RDM_mask_A_all(:,:,t) = RDM_mask_A;
    RDM_mask_B_all(:,:,t) = RDM_mask_B;

    %感觉cfar没最佳，RDM_mask_A、B存在几个点
    cfar_ranges_A = cfar_ranges_A + 1;
    cfar_ranges_B = cfar_ranges_B + 1;
    cfar_dopps_A = cfar_dopps_A + 1;
    cfar_dopps_B = cfar_dopps_B + 1;
    TrgtR = Range(cfar_dopps_A);
    TrftV = Velovity(cfar_ranges_A);

    % 获取对应目标在波束 A 和 B 中的强度
    intensity_A = abs(FFT_y1out(cfar_ranges_A, cfar_dopps_A));
    intensity_B = abs(FFT_y2out(cfar_ranges_B, cfar_dopps_B));
        
    % 计算和（Σ）和差（Δ）
    sum_val = intensity_A + intensity_B;
    diff_val = intensity_A - intensity_B;
        
    % 计算和差比 Δ/Σ
    sum_diff_ratio = diff_val / sum_val;
        
    % 根据和差比估计角度（lookup_angle_from_sumdiffratio为查找表或者说映射函数）
    TrgtAngle = lookup_angle_from_sumdiffratio(AB_ybili, sum_diff_ratio, theta, theta_min, theta_max); 
    
    %存下测距测速测角结果
    TrgtInform = [TrgtR;TrftV;TrgtAngle];
    Detect_Result(:,t) = TrgtInform;
end

% 航迹解算
r_all = Detect_Result(1,:);
theta_all = Detect_Result(3,:)+90;
xk_out = [r_all.*cosd(theta_all);r_all.*sind(theta_all)];

RMSE_R_ave = mean(abs(Detect_Result(1,:)-Zk(1,:)));
RMSE_V_ave = mean(abs(Detect_Result(2,:)-Zk(2,:)));
RMSE_phi_ave = mean(abs(Detect_Result(3,:)-Zk(3,:)));

fprintf('平均测距误差%0.2f m \n',RMSE_R_ave);
fprintf('平均测速误差%0.3f m/s \n',RMSE_V_ave);
fprintf('平均测角误差%0.4f ° \n',RMSE_phi_ave);

%% 绘图
figure; hold on;
plot(d_LinearArray,zeros(1,M),'g^','LineWidth',1.1);
plot(Xk(1,:),Xk(2,:),'b--','LineWidth',1.1);
plot(xk_out(1,:),xk_out(2,:),'rx','LineWidth',1.1);
legend('雷达位置','真实航迹','点迹估计结果','Location','southeast');
xlabel('X（m)');ylabel('Y（m)');title('航迹检测结果');

figure; hold on;axis equal;
plot(Xk(1,:),Xk(2,:),'b--','LineWidth',1.1);
plot(xk_out(1,:),xk_out(2,:),'rx','LineWidth',1.1);
legend('真实航迹','点迹估计结果','Location','southeast');
xlabel('X（m)');ylabel('Y（m)');title('航迹放大图');

figure;
mesh(Range,Velovity,abs(FFT_y1out_all(:,:,40)));hold on;
xlabel('距离（m)');ylabel('速度（m/s)');title('MTD-波束A-距离速度检测');
set(gca, 'YDir', 'reverse');

figure;
mesh(Range,Velovity,abs(FFT_y2out_all(:,:,40)));hold on;
xlabel('距离（m)');ylabel('速度（m/s)');title('MTD-波束B-距离速度检测');
set(gca, 'YDir', 'reverse');

figure;
mesh(Range,Velovity,abs(RDM_mask_A_all(:,:,40)));hold on;
xlabel('距离（m)');ylabel('速度（m/s)');title('CFAR-波束A-距离速度检测');
set(gca, 'YDir', 'reverse');

figure;
mesh(Range,Velovity,abs(RDM_mask_B_all(:,:,40)));hold on;
xlabel('距离（m)');ylabel('速度（m/s)');title('CFAR-波束B-距离速度检测');
set(gca, 'YDir', 'reverse');

plotTrajectory(Detect_Result);  % 绘制航迹的动态显示


