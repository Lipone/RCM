% 指定文件夹路径
folder_path = './dataset';

% 获取所有 .mat 文件列表
file_list = dir(fullfile(folder_path, '*.mat'));
num_files = length(file_list);

% 预定义样本数量（可从第一个文件中提前读取）
temp = load(fullfile(folder_path, file_list(1).name));
data_example = temp.data_chunk';
[num_rows, ~] = size(data_example);

% 预分配 index 矩阵
index = zeros(num_files, num_rows, num_rows);

% 并行处理
parfor i = 1:num_files
    % 构建完整路径
    file_path = fullfile(folder_path, file_list(i).name);
    
    % 加载数据
    temp = load(file_path);
    if isfield(temp, 'data_chunk')
        data = temp.data_chunk';
    else
        warning('文件 %s 中没有名为 data_chunk 的变量', file_list(i).name);
        continue;
    end

    % 初始化每个文件的临时矩阵
    index_temp = zeros(num_rows, num_rows);
    for row = 1:num_rows
        for col = 1:num_rows
            if row ~= col
                index_temp(row, col) = MCM(data(row,:), data(col,:), 5, 10);
            end
        end
    end

    % 存入主 index 矩阵
    index(i, :, :) = index_temp;
end

% 求平均矩阵
MCM_index = squeeze(mean(index, 1));  % [num_rows x num_rows]