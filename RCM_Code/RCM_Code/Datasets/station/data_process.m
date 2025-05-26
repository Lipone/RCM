%data 
clear;
clc;
% 加载原始数据
% 1. 读取CSV文件

load('temperature-v1.mat');
% 检查变量名
whos  % 查看当前工作区变量，确定主变量名称

% 创建输出文件夹
output_folder = './dataset2';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 假设主变量名为 data（你可以根据实际情况修改）
rows_per_file = 50000;
total_rows = size(dataMatrix, 1);
num_files = ceil(total_rows / rows_per_file);

for i = 1:num_files
    start_idx = (i-1)*rows_per_file + 1;
    end_idx = min(i*rows_per_file, total_rows);
    data_chunk = dataMatrix(start_idx:end_idx, :);

    % 构建保存路径和文件名
    filename = sprintf('%s/temperature_part_%03d.mat', output_folder, i);

    % 保存变量为 data（可根据需要修改变量名）
    save(filename, 'data_chunk');
end

disp('拆分完成并保存到 ./dataset2 文件夹中。');