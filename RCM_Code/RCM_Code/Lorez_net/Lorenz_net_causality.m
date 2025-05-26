% 参数设置
clc;
clear;
sigma = 10;
rho = 28; 
beta = 8/3;
num_nodes = 10;  % Lorenz 系统的数量

% 定义耦合强度和连接概率的列表
coupling_strengths = [0.1,0.5,1];

p=6;k=5;Thei=10;

connection_probs = [0.1, 0.3, 0.5];

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
                    NCMindex(run,j,i)=predict_causality(data1,data2);%NCM
                    GCindex(run,j,i)=GCmy(data1(1,1:7000),data2(1,1:7000),p);%GC
                    TEindex(run,j,i)=TEmy(data1(1,1:7000),data2(1,1:7000),p,k,Thei);%TE
                    CCMindex(run,j,i)=MCM(data1(1,1:7000),data2(1,1:7000),p,Thei);%CCM
                    end
                end
            end
        end
%         avg_auc_NCM = plot_avg_roc_with_std(predict_list, ground_truth_list, 'b', 'My Classifiers');
        
        save(sprintf("NEW_output_coupling_%d_prob_%.1f.mat", coupling_strength, connection_prob), ...
            "connectivityMatrix","NCMindex","GCindex","TEindex","CCMindex");

    end
end
% 绘制带有标准差区域的平均ROC曲线
% avg_auc = plot_avg_roc_with_std(predict_list, ground_truth_list, 'b', 'My Classifiers');
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
    fill([mean_fpr fliplr(mean_fpr)], [mean_tpr+std_tpr fliplr(mean_tpr-std_tpr)], ...
        color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    xlim([0 1]);
    ylim([0 1]);
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title(['Average ROC: ' name]);
    grid on;
    
    % 显示AUC在图上的文本
    text(0.6, 0.2, sprintf('AUC = %.3f', avg_auc), 'FontSize', 12, 'Color', color);

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
