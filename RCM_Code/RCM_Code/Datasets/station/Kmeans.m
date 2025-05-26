% 假设你有一个矩阵 A
A = data;  % 示例矩阵，可以换成你自己的

% 去除对角线元素并转换为一维
n = size(A, 1);
A_no_diag = A(~eye(n));  % 去除对角线元素，返回列向量

% 聚类为两类（使用K-means）
k = 2;  % 聚类数目
[idx, C] = kmeans(A_no_diag, k);  % idx 为分类结果，C 为中心

% 显示结果
disp('聚类标签:');
disp(idx);
disp('聚类中心:');
disp(C);

% 统计每一类的数量
count_class1 = sum(idx == 1);
count_class2 = sum(idx == 2);

% 显示统计结果
fprintf('第1类数量: %d\n', count_class1);
fprintf('第2类数量: %d\n', count_class2);