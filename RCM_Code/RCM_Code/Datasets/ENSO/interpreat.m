% 加载原始数据
load('ENSO.mat');

% 插值的目标点数量（例如插值到原数据的2倍）
interpolation_factor = 10;
num_points = interpolation_factor * size(nino, 1); 

% 原始点的索引
x_original = 1:size(nino, 1);

% 新插值点的索引（在原范围内均匀分布）
x_new = linspace(1, size(nino, 1), num_points);

% 初始化插值后的数据矩阵
interpolated_data = zeros(num_points, size(nino, 2));

% 定义插值方法列表
methods = {'linear', 'spline', 'pchip', 'nearest', 'cubic'};

% 循环遍历每种插值方法
for m = 1:length(methods)
    interp_method = methods{m};  % 获取当前插值方法
    
    % 对每一列进行插值
    for i = 1:size(nino, 2)
        interpolated_data(:, i) = interp1(x_original, nino(:, i), x_new, interp_method);
    end
    
    % 构建文件名，并保存插值后的数据
    filename = sprintf('interpolated_ENSO_%s.mat', interp_method);
    save(filename, 'interpolated_data');
    
    % 显示保存文件路径的提示
    fprintf('Interpolated data saved to %s\n', filename);
end
