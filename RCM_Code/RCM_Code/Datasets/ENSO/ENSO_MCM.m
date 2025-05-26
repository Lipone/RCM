clear;
clc;
% 加载原始数据
load('ENSO2.mat');

% 插值的目标点数量（例如插值到原数据的2倍）
interpolation_factor = {20};

methods = { 'spline'};




for t=1:length(interpolation_factor)
     
    num_points = interpolation_factor{t} * size(nino, 1); 
    
    % 原始点的索引
    x_original = 1:size(nino, 1);
    
    % 新插值点的索引（在原范围内均匀分布）
    x_new = linspace(1, size(nino, 1), num_points);
    
    % 初始化插值后的数据矩阵
    interpolated_data = zeros(num_points, size(nino, 2));
    
    % 循环遍历每种插值方法
    for k = 1:length(methods)
        interp_method = methods{k};  % 获取当前插值方法
        
        % 对每一列进行插值
        for i = 1:size(nino, 2)
            interpolated_data(:, i) = interp1(x_original, nino(:, i), x_new, interp_method);
        end 
    end
    
           
           
            data=interpolated_data(:,:)';
            [~, num_columns] = size(data);
            parfor row=1:4
                for col=1:4
                    if row ~= col
                   index(row,col)=MCM(data(row,:),data(col,:),5,10);
                    end
                end
            end

            
end

      