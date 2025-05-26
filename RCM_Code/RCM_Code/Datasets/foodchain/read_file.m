clear;
clc;
% 打开文件
fileID = fopen('output_total.txt', 'r');

% 初始化存储 all_indices 的 cell 数组
all_indices = {};
count = 1;  % 计数器

% 逐行读取文件内容
while ~feof(fileID)
    line = fgetl(fileID);  % 读取一行
    
    if contains(line, 'Combination')  % 找到组合开始的标识
        % 读取 Interpolation Factor, Scale, Method, Norm1, Norm2
        line = fgetl(fileID);  % 读取插值因子、尺度、方法、Norm1、Norm2 的行
        tokens = regexp(line, 'Interpolation Factor: (\d+), Scale: (\d+), Method: (\w+), Norm1: (\d), Norm2: (\d)', 'tokens');
        tokens = tokens{1};  % 获取解析到的 tokens
        
        % 将解析的结果存储到结构体中
        all_indices{count}.interpolation_factor = str2double(tokens{1});
        all_indices{count}.scale = str2double(tokens{2});
        all_indices{count}.interp_method = tokens{3};
        all_indices{count}.norm1 = logical(str2double(tokens{4}));  % 将 '1'/'0' 转为布尔值
        all_indices{count}.norm2 = logical(str2double(tokens{5}));  % 将 '1'/'0' 转为布尔值
        
        % 读取 "Index Matrix:" 行
        line = fgetl(fileID);  % 跳过 "Index Matrix:" 行
        
        % 读取 index 矩阵
        index_matrix = [];
        for row = 1:4  % 假设 index 矩阵是 4x4
            line = fgetl(fileID);  % 读取每一行
            row_data = str2num(line);  % 将字符串转为数值数组
            index_matrix = [index_matrix; row_data];  % 将该行添加到矩阵中
        end
        
        % 将矩阵存储到结构体中
        all_indices{count}.index = index_matrix;
        
        % 计数器加一
        count = count + 1;
        
        % 跳过分隔符行 "------------------------------------------------------"
        fgetl(fileID);
    end
end

% 关闭文件
fclose(fileID);

% 输出 all_indices 以验证结果
% disp(all_indices);

%% 

% 输出所有存储的 index 到文件
fileID = fopen('output_target.txt', 'w');  % 打开文件
for i = 1:length(all_indices)
    index_matrix = all_indices{i}.index(2:end, 2:end);
    
    % 展开 index 矩阵并找到前 4 个最大值
    flattened_index = index_matrix(:);  % 展平矩阵为向量
    [sorted_values, ~] = sort(flattened_index, 'descend');  % 从大到小排序
    top_4_values = sorted_values(1:4);  % 获取前4个最大值
    
    % 检查是否 index(2,4)、index(3,2)、index(4，1) 是矩阵中最大的4个数
    if ismember(index_matrix(1,3), top_4_values) && ...
       ismember(index_matrix(2,1), top_4_values) && ...
       ismember(index_matrix(1,2), top_4_values) && ...
       ismember(index_matrix(3,1), top_4_values)
   
        % 输出符合条件的 Combination 和矩阵到文件
        fprintf(fileID, '\nCombination %d:\n', i);
        fprintf(fileID, 'Interpolation Factor: %d, Scale: %d, Method: %s, Norm1: %s, Norm2: %s\n', ...
            all_indices{i}.interpolation_factor, all_indices{i}.scale, all_indices{i}.interp_method, ...
            num2str(all_indices{i}.norm1), num2str(all_indices{i}.norm2));
        
        % 输出 Index 矩阵
        fprintf(fileID, 'Index Matrix:\n');
        for row = 1:size(index_matrix, 1)
            fprintf(fileID, '%f ', index_matrix(row, :));  % 写入该行
            fprintf(fileID, '\n');  % 换行
        end
        disp(all_indices{i}.index);
        fprintf(fileID, '------------------------------------------------------\n');
    end
end

fclose(fileID);  % 关闭文件

%% 

% 初始化用于存储满足条件的矩阵及其差值
valid_indices = {};
differences = [];

% 输出所有存储的 index 到文件
fileID = fopen('output_target.txt', 'w');  % 打开文件
for i = 1:length(all_indices)
    index_matrix = all_indices{i}.index(2:end, 2:end);  % 去掉第一行和第一列
    
    % 展开 index 矩阵并找到前 5 个最大值
    flattened_index = index_matrix(:);  % 展平矩阵为向量
    [sorted_values, ~] = sort(flattened_index, 'descend');  % 从大到小排序
    top_4_values = sorted_values(1:4);  % 获取前4个最大值
    
    % 检查是否 index(2,4)、index(3,2)、index(4,1) 是矩阵中最大的4个数
    if ismember(index_matrix(1,3), top_4_values) && ...  % 原 index(2,4)
       ismember(index_matrix(2,1), top_4_values) && ...  % 原 index(3,2)
       ismember(index_matrix(1,2), top_4_values) && ...  % 原 index(1,2)
       ismember(index_matrix(3,1), top_4_values)         % 原 index(4,1)
   
        % 找到第四大和第五大的值
        fourth_largest = sorted_values(4);
        fifth_largest = sorted_values(5);
        first_largest=sorted_values(1);
        diff = (fourth_largest - fifth_largest)/first_largest;  % 计算差值
        
        % 存储满足条件的矩阵和差值
        valid_indices{end+1}.index_matrix = index_matrix;
        valid_indices{end}.combination_number = i;  % 记录矩阵的组合编号
        valid_indices{end}.interpolation_factor = all_indices{i}.interpolation_factor;
        valid_indices{end}.scale = all_indices{i}.scale;
        valid_indices{end}.interp_method = all_indices{i}.interp_method;
        valid_indices{end}.norm1 = all_indices{i}.norm1;
        valid_indices{end}.norm2 = all_indices{i}.norm2;
        valid_indices{end}.diff = diff;  % 记录差值
    end
end

% 假设 valid_indices 是一个 cell 数组，需要先解包
diff_values = cellfun(@(x) x.diff, valid_indices);  % 提取所有 diff 字段的值

% 根据差值进行排序
[~, sorted_idx] = sort(diff_values, 'descend');  % 根据差值进行降序排序

% 找到前3个差值最大的矩阵
top_3_indices = sorted_idx(1:min(10, length(sorted_idx)));  % 取前三个


% 输出差值最大的三个矩阵
for j = 1:length(top_3_indices)
    idx = top_3_indices(j);
    combination_number = valid_indices{idx}.combination_number;
    index_matrix = valid_indices{idx}.index_matrix;
    
    % 输出符合条件的 Combination 和矩阵到文件
    fprintf(fileID, '\nTop %d: Combination %d:\n', j, combination_number);
    fprintf(fileID, 'Interpolation Factor: %d, Scale: %d, Method: %s, Norm1: %s, Norm2: %s\n', ...
        valid_indices{idx}.interpolation_factor, valid_indices{idx}.scale, valid_indices{idx}.interp_method, ...
        num2str(valid_indices{idx}.norm1), num2str(valid_indices{idx}.norm2));
    
    fprintf(fileID, 'Index Matrix:\n');
    for row = 1:size(index_matrix, 1)
        fprintf(fileID, '%.6f ', index_matrix(row, :));  % 输出 6 位小数
        fprintf(fileID, '\n');  % 换行
    end
    
    % 输出第四大和第五大的差值
    fprintf(fileID, 'Difference between 4th and 5th largest: %.6f\n', valid_indices{idx}.diff);
    fprintf(fileID, '------------------------------------------------------\n');
end

fclose(fileID);  % 关闭文件


