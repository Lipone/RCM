% 输出每种 norm1 和 norm2 组合的前3个最佳结果
clc;
load('best_result.mat')
fields = fieldnames(best_results);
for i = 1:numel(fields)
    combination = fields{i};
    fprintf('\n%s 的前3个最佳结果:\n', combination);
    for j = 1:3
        fprintf('第 %d 名: \n', j);
        fprintf('差值: %f\n', best_results.(combination).top3_diffs(j));
        fprintf('插值因子: %d\n', best_results.(combination).top3_ts(j));
        fprintf('尺度因子: %d\n', best_results.(combination).top3_ss(j));
        disp('对应的index矩阵为:');
        disp(best_results.(combination).top3_indices{j});
    end
end