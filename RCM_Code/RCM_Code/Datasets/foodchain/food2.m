clear;
clc;
% 加载原始数据
load ("foodchain.mat");

nino= time_FCCopy(:,:);

% 插值的目标点数量（例如插值到原数据的2倍）
interpolation_factor = {5,10};
scale={5,6,10};
% 定义插值方法列表
% methods = {'linear', 'spline', 'pchip', 'nearest', 'cubic'};
methods = { 'spline'};
count = 1;
% 初始化一个结构数组来存储每种 norm1 和 norm2 组合的最佳结果
% best_results = struct();

% % 初始化最大差值及其对应的轮次信息
% max_diff = -inf;
% best_t = 0;
% best_s = 0;
% best_index = [];  % 用于存储最大差值对应的index矩阵

% 遍历 norm1 和 norm2 的 true 和 false 组合
for norm1_val = [true, false]   % norm1_val 依次取 true 和 false
    for norm2_val = [true, false]  % norm2_val 依次取 true 和 false
        
        % 跳过 norm1_val 和 norm2_val 都为 false 的组合
        % if ~norm1_val && ~norm2_val
        %     continue;  % 跳过当前循环，进入下一个组合
        % end
    
    % 初始化一个数组来存储前3个最大差值及其对应的信息
    % top3_diffs = -inf(1, 3);  % 存储前3个差值，初始值为 -inf
    % top3_ts = zeros(1, 3);    % 存储对应的插值因子
    % top3_ss = zeros(1, 3);    % 存储对应的尺度因子
    % top3_indices = cell(1, 3); % 存储对应的 index 矩阵
    % 

for t=1:length(interpolation_factor)
        
    num_points = interpolation_factor{t} * size(nino, 1); 
    
    % 原始点的索引
    x_original = 1:size(nino, 1);
    
    % 新插值点的索引（在原范围内均匀分布）
    x_new = linspace(1, size(nino, 1), num_points);
    
    % 初始化插值后的数据矩阵
    interpolated_data = zeros(num_points, size(nino, 2));
    
    % 循环遍历每种插值方法
    % for k = 1:length(methods)
    %     interp_method = methods{k};  % 获取当前插值方法
    % 
    %     % 对每一列进行插值
    %     for i = 1:size(nino, 2)
    %         interpolated_data(:, i) = interp1(x_original, nino(:, i), x_new, interp_method);
    %     end 
    % end
    
    
    % 循环遍历放缩
    for s = 1:length(scale)
        for m = 1:length(methods)
            interp_method = methods{m};  % 获取当前插值方法
            for i = 1:size(nino, 2)
                interpolated_data(:, i) = interp1(x_original, nino(:, i), x_new, interp_method);
            end 
            fprintf('interpolation_factor: %d\t\t', interpolation_factor{t});
            fprintf('scale: %d\t', scale{s});
            fprintf('interp_method: %s\n', interp_method);
            fprintf('norm1: %s\t', num2str(norm1_val));
            fprintf('norm2: %s\n', num2str(norm2_val));
            data=interpolated_data(:,:)';
            [~, num_columns] = size(data);
            for row=1:4
                for col=1:4
                    if row ~= col
                   index(row,col)=NCM(data(row,:),data(col,:),1,num_columns, scale{s}, norm1_val, norm2_val);
                    end
                end
            end

            % % 找出 index 矩阵的前 6 个最大值
            % all_values = index(:);  % 将 index 矩阵展平为向量
            % sorted_values = sort(all_values, 'descend');  % 对所有值从大到小排序
            % top_6_values = sorted_values(1:6);  % 提取前 8 个最大值
            % 
            % % 计算 index[3,4] - index[4,2]
            % if index(4, 2) <= 0.3    && ismember(index(3, 4), top_6_values)
            %     diff = index(3, 4) - index(4, 2);
            %     fprintf('index[3,4] - index[4,2] = %f\n', diff);
            % 
            % % 检查当前差值是否属于前3个最大值
            %     if diff > top3_diffs(3)  % 只比较最小的top3差值
            %         % 如果当前差值更大，更新前3个差值
            %         top3_diffs(3) = diff;
            %         top3_ts(3) = interpolation_factor{t};
            %         top3_ss(3) = scale{s};
            %         top3_indices{3} = index;
            % 
            %         % 按照差值重新排序（从大到小）
            %         [top3_diffs, sort_idx] = sort(top3_diffs, 'descend');
            %         top3_ts = top3_ts(sort_idx);
            %         top3_ss = top3_ss(sort_idx);
            %         top3_indices = top3_indices(sort_idx);
            %     end
            % 
            % end

            disp(index);
            % 打印分割线
            fprintf('------------------------------------------------------\n');
            % 存储该 index 矩阵及其对应的参数信息
            all_indices{count}.index = index;
            all_indices{count}.interpolation_factor = interpolation_factor{t};
            all_indices{count}.scale = scale{s};
            all_indices{count}.interp_method = interp_method;
            all_indices{count}.norm1 = norm1_val;
            all_indices{count}.norm2 = norm2_val;
            
            % 计数器加一
            count = count + 1;
        end
    end
end

       % 存储当前 norm1 和 norm2 组合的前3个最佳结果
        % combination_key = sprintf('norm1_%d_norm2_%d', norm1_val, norm2_val);
        % best_results.(combination_key).top3_diffs = top3_diffs;
        % best_results.(combination_key).top3_ts = top3_ts;
        % best_results.(combination_key).top3_ss = top3_ss;
        % best_results.(combination_key).top3_indices = top3_indices;
    end
end

% 输出每种 norm1 和 norm2 组合的前3个最佳结果
% fields = fieldnames(best_results);
% for i = 1:numel(fields)
%     combination = fields{i};
%     fprintf('\n%s 的前3个最佳结果:\n', combination);
%     for j = 1:3
%         fprintf('第 %d 名: \n', j);
%         fprintf('差值: %f\n', best_results.(combination).top3_diffs(j));
%         fprintf('插值因子: %d\n', best_results.(combination).top3_ts(j));
%         fprintf('尺度因子: %d\n', best_results.(combination).top3_ss(j));
%         disp('对应的index矩阵为:');
%         disp(best_results.(combination).top3_indices{j});
%     end
% end

% Step 1: 读取文件并查找现有的最大 Combination 编号
fileID = fopen('output_total.txt', 'r');
max_combination = 0;  % 初始化为 0

if fileID ~= -1  % 检查文件是否存在
    while ~feof(fileID)
        line = fgetl(fileID);  % 读取一行
        if contains(line, 'Combination')
            tokens = regexp(line, 'Combination (\d+):', 'tokens');  % 提取 Combination 编号
            if ~isempty(tokens)
                current_combination = str2double(tokens{1}{1});  % 转换为数值
                max_combination = max(max_combination, current_combination);  % 更新最大值
            end
        end
    end
    fclose(fileID);  % 关闭文件
else
    disp('File does not exist, starting from Combination 1.');
end
% 输出所有存储的 index
fileID = fopen('output_total.txt', 'a');
for i = 1:length(all_indices)
    fprintf(fileID, '\nCombination %d:\n', max_combination + i);
    fprintf(fileID, 'Interpolation Factor: %d, Scale: %d, Method: %s, Norm1: %s, Norm2: %s\n', ...
        all_indices{i}.interpolation_factor, all_indices{i}.scale, all_indices{i}.interp_method, ...
        num2str(all_indices{i}.norm1), num2str(all_indices{i}.norm2));
    fprintf(fileID, 'Index Matrix:\n');
    index_matrix = all_indices{i}.index;
    
    % 逐行写入 index 矩阵
    for row = 1:size(index_matrix, 1)
        fprintf(fileID, '%f ', index_matrix(row, :));  % 写入该行
        fprintf(fileID, '\n');  % 换行
    end
    
    fprintf(fileID, '------------------------------------------------------\n');

end
fclose(fileID);
