%%
%测试训练长度
clear;
clc;
CList = [153,153,253; 255,153,154]./255;

adj_matrix=[ 0,0,0,0,1;
    1,0,0,0,0;
    0,1,0,0,0;
    0,0,1,0,0;
    0,0,0,1,0;
    ];
sigma = 10;
rho = 28;
beta = 8/3;
numNodes=5;
coupling_strengths = [0.1];
tspan = linspace(0, 1000, 10000);
trails=10;
NCMindex=zeros(trails,numNodes,numNodes);
outter_noises=[0];
inner_noises=[0];
resSize=[10,20,30,50,100,200,300,400,500,1000];
result1 = zeros(trails,1);  % 预分配
result2 = zeros(trails,1);
dim=3;

all_ones=cell(length(resSize), 1);
all_zeros=cell(length(resSize), 1);
%%
step=1;
for res=resSize
    parfor trail = 1:trails
        dim2=res;
        sigma=0.01;
        Win1 = -1 + 2*rand(1000,dim2);
        adj1 = zeros(1000,dim2);
        for m=1:1000
            for n=1:dim2
                if(rand(1,1)<sigma)  
                    adj1(m,n)=1;  
                end
            end
        end
        Win_predict = adj1.*Win1;

        filename = sprintf('ESN_weights_EmbSize_200_0%d.mat', trail);  % 格式化文件名
        filepath = fullfile('Weights', filename);           % 拼接文件路径
        S = load(filepath, 'Win', 'Wres');
        Win1 = S.Win;
        Wres1 = S.Wres;

        filename = sprintf('ESN_weights_EmbSize_%d_0%d.mat', res,trail);  % 格式化文件名
        filepath = fullfile('Weights', filename);           % 拼接文件路径
        S = load(filepath, 'Win', 'Wres');
        Win2= S.Win;
        Wres2 = S.Wres;                    % 加载权重
%         Win_predict=S.Win_predict;


        % 初始条件 (10 个 Lorenz 系统，每个系统 3 个变量)
        x0 = rand(3 * numNodes, 1);  % 生成初始条件，长度为 3*numNodes
        % 计算系统动态
        [t, X] = ode45(@(t,x) lorenz_system(t, x, 10, rho, beta, coupling_strengths, numNodes, adj_matrix,0), tspan, x0);
        X = X + normrnd(0, 0, size(X));
        % 初始化返回矩阵
        num_time_steps = length(t);  % 时间步数
        node_states = zeros(numNodes, num_time_steps, 3);  % 三维矩阵，节点数 x 时间步数 x 状态数

        % 填充返回矩阵，每个节点有 x, y, z 三个状态
        for i = 1:numNodes
            node_states(i, :, 1) = X(:, 3*i-2);  % x 状态
            node_states(i, :, 2) = X(:, 3*i-1);  % y 状态
            node_states(i, :, 3) = X(:, 3*i);    % z 状态
        end

        data1=squeeze(node_states(1,:,:))';
        data2=squeeze(node_states(2,:,:))';
        % 执行计算
        r1 = RCM(data1,data2,3,res,10000,5,Win2,Win1,Wres2,true,false,Wres1,Win_predict,false,trail);
        r2 = RCM(data2,data1,3,res,10000,5,Win2,Win1,Wres2,true,false,Wres1,Win_predict,false,trail);

        % 存入临时变量（每个 worker 写自己，不冲突）
        result1(trail) = r1;
        result2(trail) = r2;
    end
    all_ones{step}=result1;
    all_zeros{step}= result2;
    step=step+1
end


%% 

%  plot lines
color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.905882352941177,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.305882352941177,0.541176470588235,0.317647058823529;
    1.0,0.5,0.0;
    0.5,0.7,0.4;
    0, 0, 0.75;];

% 初始化图像
figure('color',[1 1 1]);
L=length(resSize)
% 第一组数据：NCM_values_for_ones
MEAN=zeros(1,L);
S=zeros(1,L);
for i = 1:L
    MEAN(:,i)=mean(all_ones{i});
    S(:,i)=std(all_ones{i});
end

x=1:1:L;
yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_ones = color_all(6,:);
color_between_ones = color_all(6,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[X→Y]', 'color', color_mean_ones);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between_ones, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;
% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_zeros)
    y = [y, all_zeros{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, color_all(8,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;

% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_ones)
    y = [y, all_ones{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, color_all(6,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;


% 第二组数据：NCM_values_for_zeros
MEAN=zeros(1,L);
S=zeros(1,L);
for i = 1:L
    MEAN(:,i)=mean(all_zeros{i});
    S(:,i)=std(all_zeros{i});
end

yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_zeros = color_all(8,:);
color_between_zeros = color_all(8,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[Y→X]', 'color', color_mean_zeros);
fill([x fliplr(x)], [yu fliplr(yl)], color_between_zeros, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;

% 设置坐标轴属性
ax = gca;
ax.LineWidth = 1.3;
ax.XAxis.Color = [0, 0, 0]; 
ax.YAxis.Color = [0, 0, 0]; 
xlabel('$\mathbf{\beta_{xy}}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');

% 设置y轴范围
% ylim([0, 0.6]);  % 将y轴范围设置为[-1, 1]，可以根据实际数据调整

% 添加图例
legend({'CCM[X→Y]', 'CCM[Y→X]'}, 'Location', 'southeast', 'TextColor', 'k', 'Box', 'off');  
legend('show');  

% 隐藏顶部和右侧坐标轴
ax.Box = 'off';
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';

% 定义耦合 Lorenz 系统的微分方程
function dxdt = lorenz_system(t, x, sigma, rho, beta, coupling_strength, numNodes, adj_matrix, noise)
if nargin < 9
    noise = 0;
end
dxdt = zeros(3 * numNodes, 1);  % 初始化导数向量

% 遍历每个节点，计算其动态方程
for i = 1:numNodes
    x_i = x(3*i-2);
    y_i = x(3*i-1);
    z_i = x(3*i);

    % Lorenz 系统方程
    dx = sigma * (y_i - x_i);
    dy = x_i * (rho - z_i) - y_i;
    dz = x_i * y_i - beta * z_i;

    % 耦合项，计算与连接的节点的耦合（根据邻接矩阵 adj_matrix）
    coupling_x = 0;
    for j = 1:numNodes
        if adj_matrix(i, j) == 1  % 如果节点 i 和节点 j 之间有连接
            x_j = x(3*j-2);
            coupling_x = coupling_x + coupling_strength * x_j;
        end
    end

    % 添加噪声项（均值为0，方差为noise的正态分布噪声）
    noise_x = normrnd(0, noise);  % x 方向的噪声
    noise_y = normrnd(0, noise);  % y 方向的噪声
    noise_z = normrnd(0, noise);  % z 方向的噪声

    % 更新当前节点的方程，添加耦合项和噪声
    dxdt(3*i-2) = dx + coupling_x + noise_x;  % x 方向加上耦合项和噪声
    dxdt(3*i-1) = dy ;  % y 方向加上噪声
    dxdt(3*i) = dz ;    % z 方向加上噪声
end
end


