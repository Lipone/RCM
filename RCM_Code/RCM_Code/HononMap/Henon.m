clear;
clc;
numNodes=8;
alpha=0.1;
sigma=0.02;

T = 1500;
x = simulate_coupled_map(numNodes, T, alpha, sigma)
%%
numNodes=8;
alpha=0.1;
sigma=0.02;
T = 1500;
trails=1;

NCMindex=zeros(trails,numNodes,numNodes);

A = zeros(numNodes);
p=5;k=7;Thei=10;
for i = 2:numNodes-1
    A(i, i-1) = 1;
    A(i, i+1) = 1;
end
ones_positions = find(A == 1 & ~eye(size(A)));  % adj_matrix 中为 1 且不在对角线的元素
zeros_positions = find(A == 0 & ~eye(size(A)));  % adj_matrix 中为 0 且不在对角线的元素

filename = sprintf('ESN_weights_01.mat');  % 格式化文件名
filepath = fullfile('Weights', filename);           % 拼接文件路径
load(filepath, 'Win', 'Wres');                      % 加载权重
for trail=1:trails
    for i=1:numNodes
        for j=1:numNodes
            if i ~= j
%                                 NCMindex(trail,j,i)=RCM(x(i,:),x(j,:),1,T,5,true,true,Win,Wres);%NCM,NCMindex(trail,j,i)表示节点i对j的因果
                NCMindex(trail,j,i)=MCM(x(i,:),x(j,:),p,Thei);%NCM,NCMindex(trail,j,i)表示节点i对j的因果
            end
        end
    end
end


for i=1:trails
    data_all=squeeze(NCMindex(i,:,:))
    if i ==1
        NCM_values_for_ones=data_all(ones_positions);
        NCM_values_for_zeros=data_all(zeros_positions);
    else
        NCM_values_for_ones=[NCM_values_for_ones,data_all(ones_positions)];  % 对应于 adj_matrix 中为 1 的 NCMindex 值
        NCM_values_for_zeros= [NCM_values_for_zeros,data_all(zeros_positions)];  % 对应于 adj_matrix 中为 0 的 NCMindex 值
    end
end
for trial =trails

    NCMindex1=squeeze(NCMindex(trial,:,:));
    B=[];
    [m,n]=size(A);
    for i=1:m
        B=[B;A(i,[1:i-1 i+1:n])];
    end
    connectivityMatrix1=B;

    A = NCMindex1;
    B=[];
    [m,n]=size(A);
    for i=1:m
        B=[B;A(i,[1:i-1 i+1:n])];
    end
    NCMindex1=B;
    label=connectivityMatrix1(:)';
    p=NCMindex1(:)';
    %%
    fig = figure;
    [NCMauc,Threshold]=plot_roc(p,label,[0.3, 0.3, 0.9],'NCM')
    % Get the AUC value as a string
    NCMauc = sprintf('%.4f', NCMauc );
    count1 = sum(NCM_values_for_ones <Threshold)
    count2 = sum(NCM_values_for_zeros >= Threshold)
    erros_for_NCM=count1+count2

    % 添加图例
    legend({['MCM (AUC=' NCMauc ')']},'Location', 'southeast', 'TextColor', 'k', 'Box', 'off');  % 'Box' 属性设置为 'off' 取消图例的边框
    legend('show');  % 显示图例
    exportgraphics(fig, ['MCM', '.jpg'], 'Resolution', 1000);
end
%%

function [auc, best_threshold] = plot_roc(predict, ground_truth, color, name)
% ROC 曲线绘制，并计算 AUC 与最佳阈值
% 输入：
%   predict: 预测概率 (列向量)
%   ground_truth: 真实标签（0 或 1）
%   color: 绘图颜色
%   name: 曲线名称（可选）

% 所有可能的阈值
%初始点为（1.0, 1.0）
x = 1.0;
y = 1.0;
%计算出ground_truth中正样本的数目pos_num和负样本的数目neg_num
pos_num = sum(ground_truth==1);
neg_num = sum(ground_truth==0);
%根据该数目可以计算出沿x轴或者y轴的步长
x_step = 1.0/neg_num;
y_step = 1.0/pos_num;
%首先对predict中的分类器输出值按照从小到大排列
[predict,index] = sort(predict);
ground_truth = ground_truth(index);
%对predict中的每个样本分别判断他们是FP或者是TP
%遍历ground_truth的元素，
%若ground_truth[i]=1,则TP减少了1，往y轴方向下降y_step
%若ground_truth[i]=0,则FP减少了1，往x轴方向下降x_step
X(1)=1;
Y(1)=1;
for i=1:length(ground_truth)
    if ground_truth(i) == 1
        y = y - y_step;
    else
        x = x - x_step;
    end
    X(i+1)=x;
    Y(i+1)=y;
end
%画出图像


plot(X, Y, 'LineWidth', 2, 'Color', color);
thresholds = unique(predict);
n = length(thresholds);

TPR = zeros(n, 1);
FPR = zeros(n, 1);
J_stat = zeros(n, 1);
F1=zeros(n, 1);

pos_num = sum(ground_truth == 1);
neg_num = sum(ground_truth == 0);

for i = 1:n
    th = thresholds(i);
    pred_binary = predict >= th;

    TP = sum((pred_binary == 1) & (ground_truth == 1));
    FP = sum((pred_binary == 1) & (ground_truth == 0));
    FN = sum((pred_binary == 0) & (ground_truth == 1));
    TN = sum((pred_binary == 0) & (ground_truth == 0));

    TPR(i) = TP / (TP + FN);  % Sensitivity
    FPR(i) = FP / (FP + TN);  % 1 - Specificity
    Precision = TP / (TP + FP + eps);  % 加 eps 避免除 0
    Recall = TP / (TP + FN + eps);
    F1(i) = 2 * Precision * Recall / (Precision + Recall + eps);
    J_stat(i) = TPR(i) - FPR(i);
end

% 计算 AUC（按 FPR 升序排序）
[FPR_sorted, idx] = sort(FPR);
TPR_sorted = TPR(idx);
auc = trapz(FPR_sorted, TPR_sorted);  % 正向积分

% 找到最大 Youden's J statistic 对应的阈值
[~, best_idx] = max(F1);
best_threshold = thresholds(best_idx);
% 绘制 ROC 曲线
hold on;
%     plot(FPR_sorted, TPR_sorted, 'LineWidth', 2, 'Color', color, 'LineStyle', '-', 'DisplayName', name);
plot(FPR(best_idx), TPR(best_idx), 'ko', 'MarkerSize', 8, 'LineWidth', 2);  % 最佳点标记
hold off;
ax = gca;
ax.LineWidth = 1.3;
ax.Box = 'off';
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';
ax.FontName = 'Times New Roman';
ax.FontSize = 20;
ax.FontWeight = 'bold';      % 设置字体加粗


xlabel('False Positive Rate', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('True Positive Rate', 'FontSize', 20, 'FontWeight', 'bold');

end


%%

% figure;
% plot(sort(NCM_values_for_ones));
% hold on;
%
% min(NCM_values_for_ones)
function x = simulate_coupled_map(n, T, alpha, sigma)
% 模拟耦合非线性动力学系统
%
% 输入：
%   n     - 节点个数
%   T     - 时间步长
%   alpha - 耦合强度（0 ~ 1）
%   sigma - 噪声强度
%
% 输出：
%   x     - 大小为 n x T 的状态矩阵，每行对应一个节点

% 初始化状态矩阵和噪声
x = zeros(n, T);
x(:, 1:2) = randn(n, 2);         % 前两步随机初始化
epsilon = sigma * randn(n, T);  % 高斯白噪声

% 动力学演化
for t = 3:T
    for i = 1:n
        if i == 1 || i == n
            % 边界节点：无耦合项
            x(i, t) = 1 - 1.4 *x(i, t-1)^2 + 0.2 * x(i, t-2) + epsilon(i, t);
        else
            % 内部节点：包含耦合
            neighbor_avg = 0.5 * alpha * (x(i-1, t-1) + x(i+1, t-1));
            self_term = (1 - alpha) * x(i, t-1);
            x(i, t) = 1 - 1.4 *(neighbor_avg + self_term)^2 + 0.2 * x(i, t-2) + epsilon(i, t);
%             x(i, t) = 1 - 1.4 *x(i, t-1)^2 + 0.2 * x(i, t-2) + epsilon(i, t)+alpha * (2/3*x(i-1, t-1)^2-1/5*sin(x(i-1, t-1)) + 2/3*x(i+1, t-1)-1/5*sin(x(i+1, t-1)));
        end
    end
end
end
