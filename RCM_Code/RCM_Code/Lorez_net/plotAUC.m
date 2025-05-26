clear;
clc;
color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.905882352941177,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.305882352941177,0.541176470588235,0.317647058823529;
    1.0,0.5,0.0;
    0.5,0.7,0.4;
    0, 0, 0.75;];
CList = [0, 0, 0.75; 1.0,0.5,0.0;[153,153,253; 255,153,154]./255];
colors= [0.3, 0.3, 0.9;0.5,0.7,0.2;0.9,0.7,0.2;0.1,0.7,0.9];

%
% load("NEW_output_coupling_1_prob_0.3.mat");
filename = '3D_CCMindex_coupling_1_prob_0.3.mat';  % 你的 .mat 文件名

% 获取不带扩展名的文件名
[~, name, ~] = fileparts(filename);
load(filename);
trails=length(CCMindex(:,1,1))
NCMindex = mat2cell(reshape(NCMindex, trails, 100), ones(1, trails), 100);
TEindex = mat2cell(reshape(TEindex, trails, 100), ones(1, trails), 100);
GCindex = mat2cell(reshape(GCindex, trails, 100), ones(1, trails), 100);
CCMindex = mat2cell(reshape(CCMindex, trails, 100), ones(1, trails), 100);
for i = 1:trails
    connectivityMatrix{i} = reshape(connectivityMatrix{i}, [100, 1]);
end
fig = figure;
% 绘制带有标准差区域的平均ROC曲线
NCMauc = plot_avg_roc_with_std(NCMindex, connectivityMatrix, colors(1,:), 'My Classifiers');
CCMauc = plot_avg_roc_with_std(CCMindex, connectivityMatrix, colors(4,:), 'My Classifiers');
TEauc = plot_avg_roc_with_std(TEindex, connectivityMatrix, colors(2,:), 'My Classifiers');
GCauc = plot_avg_roc_with_std(GCindex, connectivityMatrix, colors(3,:), 'My Classifiers');
% Get the AUC value as a string
NCMauc = sprintf('%.4f', NCMauc );
CCMauc = sprintf('%.4f', CCMauc );
GCauc = sprintf('%.4f', GCauc );
TEauc = sprintf('%.4f', TEauc );
% 设置 X 轴和 Y 轴范围及刻度
xlim([0 1]);
xticks(0:0.2:1);
ylim([0 1]);
yticks(0:0.2:1);

% 获取当前坐标轴对象
ax = gca;

% 设置只显示 X 轴和 Y 轴
ax.XAxis.Visible = 'on';       % 显示 X 轴
ax.YAxis.Visible = 'on';       % 显示 Y 轴
ax.Box = 'off';                % 隐藏顶部和右侧框线

% 设置字体属性
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';        % 字体加粗
ax.FontSize = 20;              % 字体大小
ax.FontName = 'Times New Roman'; % 设置字体为 Times New Roman

% 添加坐标轴标签
xlabel('False Positive Rate', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('True Positive Rate', 'FontSize', 20, 'FontWeight', 'bold');


% 添加图例
legend({['RCM (AUC=' NCMauc ')'],['MCM (AUC=' CCMauc ')'],['TE (AUC=' TEauc ')'],['GC (AUC=' GCauc ')']},'Location', 'southeast', 'TextColor', 'k', 'Box', 'off','FontWeight', 'bold');  % 'Box' 属性设置为 'off' 取消图例的边框
legend('show');  % 显示图例


exportgraphics(fig, [name, '.jpg'], 'Resolution', 1000);


%%

% load("output_coupling_1_prob_0.3.mat")
filename = 'NEW_output_coupling_1.000000e-01_prob_0.1.mat';  % 你的 .mat 文件名

% 获取不带扩展名的文件名
[~, name, ~] = fileparts(filename);
load(filename);
for dim =1:5
    figure;
    % load("output_coupling_5.000000e-01_prob_0.3.mat");
    % load NEW_output_coupling_1_prob_0.5.mat
    NCMindex1=squeeze(NCMindex(dim,:,:));
    TEindex1=squeeze(TEindex(dim,:,:));
    GCindex1=squeeze(GCindex(dim,:,:));
    CCMindex1=squeeze(CCMindex(dim,:,:));
    connectivityMatrix1=connectivityMatrix{dim}
    A = connectivityMatrix1;
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

    A = GCindex1;
    B=[];
    [m,n]=size(A);
    for i=1:m
        B=[B;A(i,[1:i-1 i+1:n])];
    end
    GCindex1=B;

    A = TEindex1;
    B=[];
    [m,n]=size(A);
    for i=1:m
        B=[B;A(i,[1:i-1 i+1:n])];
    end
    TEindex1=B;

    A = CCMindex1;
    B=[];
    [m,n]=size(A);
    for i=1:m
        B=[B;A(i,[1:i-1 i+1:n])];
    end
    CCMindex1=B;

    label=connectivityMatrix1(:)';
    p=NCMindex1(:)';
    NCMauc=plot_roc(p,label,[0.3, 0.3, 0.9],'NCM')
    % NCMaucAAA = plot_avg_roc_with_std(NCMindex, connectivityMatrix, colors(1,:), 'My Classifiers')
    hold on

    p=CCMindex1(:)';
    CCMauc=plot_roc(p,label,[0.5,0.7,0.2],'CCM')
    hold on

    p=GCindex1(:)';
    GCauc=plot_roc(p,label,[0.9,0.7,0.2],'GC')
    hold on

    p=TEindex1(:)';
    TEauc=plot_roc(p,label,[0.1,0.7,0.9],'TE')
    hold on

    % Get the AUC value as a string
    NCMauc = sprintf('%.4f', NCMauc );
    CCMauc = sprintf('%.4f', CCMauc );
    GCauc = sprintf('%.4f', GCauc );
    TEauc = sprintf('%.4f', TEauc );

    % 添加图例
    legend({['RCM (AUC=' NCMauc ')'],['MCM (AUC=' CCMauc ')'],['GC (AUC=' GCauc ')'],['TE (AUC=' TEauc ')']},'Location', 'southeast', 'TextColor', 'k', 'Box', 'off');  % 'Box' 属性设置为 'off' 取消图例的边框
    legend('show');  % 显示图例
end
function  auc = plot_roc( predict, ground_truth ,color,name )
% INPUTS
%  predict       - 分类器对测试集的分类结果
%  ground_truth - 测试集的正确标签,这里只考虑二分类，即0和1
% OUTPUTS
%  auc            - 返回ROC曲线的曲线下的面积

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
% plot(X,Y,'-ro','LineWidth',2,'MarkerSize',3);

plot(X, Y, 'LineWidth', 2, 'Color', color, 'LineStyle', '--');

% 设置 X 轴和 Y 轴范围及刻度
xlim([0 1]);
xticks(0:0.2:1);
ylim([0 1]);
yticks(0:0.2:1);

% 获取当前坐标轴对象
ax = gca;

% 设置只显示 X 轴和 Y 轴
ax.XAxis.Visible = 'on';       % 显示 X 轴
ax.YAxis.Visible = 'on';       % 显示 Y 轴
ax.Box = 'off';                % 隐藏顶部和右侧框线

% 设置字体属性
ax.FontWeight = 'bold';        % 字体加粗
ax.FontSize = 15;              % 字体大小
ax.FontName = 'Times New Roman'; % 设置字体为 Times New Roman

% 添加坐标轴标签
xlabel('False Positive Rate', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('True Positive Rate', 'FontSize', 15, 'FontWeight', 'bold');


%计算小矩形的面积,返回auc
auc = -trapz(X,Y);
end
%%

function avg_auc = plot_avg_roc_with_std(predict_list, ground_truth_list, color, name)
% INPUTS
% predict_list        - 包含多个分类器的预测值的列表，大小为 [num_classifiers x num_samples]
% ground_truth_list   - 包含多个分类器的真实标签的列表，大小为 [num_classifiers x num_samples]
% color               - ROC曲线的颜色
% name                - 曲线的名称

% 定义统一的假阳性率范围
num_points = 1000;
mean_fpr = linspace(0, 1, num_points); % 定义从0到1的均匀分布的假阳性率范围
tprs = zeros(length(predict_list), num_points); % 存储每个分类器的插值TPR
aucs = zeros(1, length(predict_list)); % 存储每个分类器的AUC

% 遍历每个分类器或交叉验证折叠
for i = 1:length(predict_list)
    predict = predict_list{i}; % 获取第 i 个分类器的预测值
    ground_truth = ground_truth_list{i}; % 获取第 i 个分类器的真实标签

    % 计算该分类器的ROC曲线
    [fpr, tpr, ~, auc] = perfcurve(ground_truth, predict, 1);
    aucs(i) = auc; % 存储AUC值
    % 去除重复的FPR值，确保插值可以顺利进行
    [fpr, unique_idx] = unique(fpr);  % 获取唯一的fpr和相应的索引
    tpr = tpr(unique_idx);  % 保留与唯一fpr对应的tpr
    % 插值到统一的假阳性率点
    tprs(i, :) = interp1(fpr, tpr, mean_fpr, 'linear', 'extrap');
end

% 计算平均TPR
mean_tpr = mean(tprs, 1); % 计算每个假阳性率下TPR的平均值
std_tpr = std(tprs, 0, 1); % 计算每个假阳性率下TPR的标准差
avg_auc = mean(aucs); % 计算AUC的平均值

% 绘制平均ROC曲线
plot(mean_fpr, mean_tpr, 'LineWidth', 2, 'Color', color);
hold on;

% 绘制标准差区域（用阴影表示）
%     fill([mean_fpr fliplr(mean_fpr)], [mean_tpr+std_tpr fliplr(mean_tpr-std_tpr)], ...
%         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlim([0 1]);
ylim([0 1]);
xlabel('False Positive Rate');
ylabel('True Positive Rate');

xlim([0 1]);
xticks(0:0.2:1);
ylim([0 1]);
yticks(0:0.2:1);
% 获取当前坐标轴对象
ax = gca;
% 隐藏顶部和右侧坐标轴
ax.Box = 'off';

% 只显示 X 轴和 Y 轴
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';
set(gca,'FontSize',20,'FontName','Times New Roman');
end