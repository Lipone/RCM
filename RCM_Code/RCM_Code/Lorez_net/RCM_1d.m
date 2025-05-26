% 参数设置
clc;
clear;
sigma = 10;
rho = 28; 
beta = 8/3;
num_nodes = 10;  % Lorenz 系统的数量

% 定义耦合强度和连接概率的列表
coupling_strengths = [0.1];

p=6;k=5;Thei=10;

connection_probs = [0.1,0.3];

% 初始条件 (10 个 Lorenz 系统，每个系统 3 个变量)
x0 = rand(3 * num_nodes, 1);  % 生成初始条件，长度为 3*num_nodes

% 时间范围
tspan = linspace(0, 500, 50000);

% 嵌套循环计算不同的耦合强度和连接概率
for coupling_strength = coupling_strengths
    for connection_prob = connection_probs
        numNodes = num_nodes;
        NCMindex=zeros(5,numNodes,numNodes);
        GCindex=zeros(5,numNodes,numNodes);
        TEindex=zeros(5,numNodes,numNodes);
        CCMindex=zeros(5,numNodes,numNodes);
        connectivityMatrix = cell(5, 1);  % 使用 cell 数组来存储每次运行的邻接矩阵
        for run= 1:5
            % 生成随机的邻接矩阵
            adj_matrix = rand(num_nodes) < connection_prob;  
            for i = 1:num_nodes
                adj_matrix(i, i) = 0;
            end
    
            % 计算系统动态
            [t, X] = ode45(@(t,x) lorenz_system(t, x, sigma, rho, beta, coupling_strength, num_nodes, adj_matrix), tspan, x0);
    
            % 初始化返回矩阵
            num_time_steps = length(t);  % 时间步数
            node_states = zeros(num_nodes, num_time_steps, 3);  % 三维矩阵，节点数 x 时间步数 x 状态数
    
            % 填充返回矩阵，每个节点有 x, y, z 三个状态
            for i = 1:num_nodes
                node_states(i, :, 1) = X(:, 3*i-2);  % x 状态
                node_states(i, :, 2) = X(:, 3*i-1);  % y 状态
                node_states(i, :, 3) = X(:, 3*i);    % z 状态
            end
    
            % 保存文件，文件名根据耦合强度和连接概率生成
            output = node_states;
            
            connectivityMatrix{run} = adj_matrix;

            
            p=6;k=5;Thei=10;
            
            for i=1:numNodes
                for j=1:numNodes
                    if i ~= j
                    data1=squeeze(output(i,:,:))';  
                    data2=squeeze(output(j,:,:))';
                    NCMindex(run,j,i)=predict_causality2(data1(1,1:7000),data2(1,1:7000));%NCM
                    GCindex(run,j,i)=GCmy(data1(1,1:7000),data2(1,1:7000),p);%GC
                    TEindex(run,j,i)=TEmy(data1(1,1:7000),data2(1,1:7000),p,k,Thei);%TE
                    CCMindex(run,j,i)=MCM(data1(:,1:7000),data2(:,1:7000),p,Thei);%CCM
                    end
                end
            end
        end
%         avg_auc_NCM = plot_avg_roc_with_std(predict_list, ground_truth_list, 'b', 'My Classifiers');
        
        save(sprintf("1D_NCMindex_coupling_%d_prob_%.1f.mat", coupling_strength, connection_prob), ...
            "connectivityMatrix","NCMindex");

    end
end
%% 
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

trails=length(NCMindex(:,1,1))
NCMindex = mat2cell(reshape(NCMindex, trails, 100), ones(1, trails), 100);

for i = 1:trails
    connectivityMatrix{i} = reshape(connectivityMatrix{i}, [100, 1]);
end
fig = figure;
% 绘制带有标准差区域的平均ROC曲线
NCMauc = plot_avg_roc_with_std(NCMindex, connectivityMatrix, colors(1,:), 'My Classifiers');

% Get the AUC value as a string
NCMauc = sprintf('%.4f', NCMauc );

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
ax.FontSize = 20;              % 字体大小
ax.FontName = 'Times New Roman'; % 设置字体为 Times New Roman

% 添加坐标轴标签
xlabel('False Positive Rate', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('True Positive Rate', 'FontSize', 20, 'FontWeight', 'bold');


% 添加图例
legend({['RCM (AUC=' NCMauc ')']},'Location', 'southeast', 'TextColor', 'k', 'Box', 'off','FontWeight', 'bold');  % 'Box' 属性设置为 'off' 取消图例的边框
legend('show');  % 显示图例

%% 

for dim =1:1
    figure;

NCMindex=squeeze(NCMindex(dim,:,:));
% CCMindex=squeeze(CCMindex(dim,:,:));
connectivityMatrix=connectivityMatrix{dim}
A = connectivityMatrix;
B=[];
[m,n]=size(A);
for i=1:m
    B=[B;A(i,[1:i-1 i+1:n])];
end
connectivityMatrix=B;

A = NCMindex;
B=[];
[m,n]=size(A);
for i=1:m
    B=[B;A(i,[1:i-1 i+1:n])];
end
NCMindex=B;

label=connectivityMatrix(:)';
p=NCMindex(:)';
NCMauc=plot_roc(p,label,[0.3, 0.3, 0.9],'NCM')
hold on

% Get the AUC value as a string
NCMauc = sprintf('%.4f', NCMauc );
% CCMauc = sprintf('%.4f', CCMauc );

% 添加图例
legend({['RCM (AUC=' NCMauc ')']},'Location', 'southeast', 'TextColor', 'k', 'Box', 'off');  % 'Box' 属性设置为 'off' 取消图例的边框
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










% 定义耦合 Lorenz 系统的微分方程
function dxdt = lorenz_system(t, x, sigma, rho, beta, coupling_strength, num_nodes, adj_matrix)
    dxdt = zeros(3 * num_nodes, 1);  % 初始化导数向量

    % 遍历每个节点，计算其动态方程
    for i = 1:num_nodes
        x_i = x(3*i-2);
        y_i = x(3*i-1);
        z_i = x(3*i);

        % Lorenz 系统方程
        dx = sigma * (y_i - x_i);
        dy = x_i * (rho - z_i) - y_i;
        dz = x_i * y_i - beta * z_i;

        % 耦合项，计算与连接的节点的耦合（根据邻接矩阵 adj_matrix）
        coupling_x = 0;
        for j = 1:num_nodes
            if adj_matrix(i, j) == 1  % 如果节点 i 和节点 j 之间有连接
                x_j = x(3*j-2);
                coupling_x = coupling_x + coupling_strength * x_j;
            end
        end

        % 更新当前节点的方程
        dxdt(3*i-2) = dx + coupling_x;  % x 方向加上耦合项
        dxdt(3*i-1) = dy;  % y 方向不耦合
        dxdt(3*i) = dz;    % z 方向不耦合
    end
end

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
