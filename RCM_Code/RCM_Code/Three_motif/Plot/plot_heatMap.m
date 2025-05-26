%% generate data
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');

%% 
load('colorData.mat')

%REAL
%fan-in
fan_in =zeros(3, 3);
fan_in (1,2)=1;
fan_in (3,2)=1;
%fan-out
fan_out =zeros(3, 3);
fan_out (1,3)=1;
fan_out (1,2)=1;
%cascade
cascade =zeros(3,3);
cascade (1,3)=1;
cascade (3,2)=1;


load("CCM_couple_strengrh_0.1.mat")
load("RCM_couple_strengrh_0.1.mat")
load("CCM_couple_strengrh_0.4.mat")
load("RCM_couple_strengrh_0.4(2).mat")



fig1 = figure;
fig1.Position = [100, 100, 800, 500];  % 设置窗口大小（像素）
t =tiledlayout(3, 5, 'Padding', 'compact', 'TileSpacing', 'compact');  % 更紧凑布局

% 依次绘制 15 张图
mat_list = {
    cascade, Cascade_RCM1, Cascade_RCM4, Cascade_CCM1, Cascade_CCM4,
    fan_out, fan_out_RCM1, fan_out_RCM4, fan_out_CCM1, fan_out_CCM4,
    fan_in, fan_in_RCM1, fan_in_RCM4, fan_in_CCM1, fan_in_CCM4
};
mat_list = mat_list';  %

for i = 1:15
    nexttile(i);
    imagesc(mat_list{i}); % 显示图像
    xticks([]); yticks([]);
    caxis([0 1]);
    axis square;
    set(gca, 'FontName', 'Times New Roman');  % 设置字体

    % 获取当前图像大小
    [rows, cols] = size(mat_list{i});
    
    % 绘制对角线
    hold on;
    plot([0.1, cols+0.9], [0.1, rows+0.9], '--k', 'LineWidth', 1,'Color', [0.5 0.5 0.5 0.3]);  % 灰色，透明度 0.5);  % 虚线对角线，黑色
    hold off;
end

colormap(othercolor('PuBu8'));
cb = colorbar;
% cb.Layout.Tile = 'east';         % ← 这是关键，放在 tiledlayout 外侧
cb.FontSize = 10;
cb.FontWeight = 'bold';
caxis([0 1]);
cb.Position = [0.96,0.1,0.012,0.82];  % [left, bottom, width, height]
cb.Layout.Tile = 'east'; 
exportgraphics(fig1, 'fig1.jpg', 'Resolution', 1500);



fig2 = figure;
fig2.Position = [100, 100, 500, 700];  % 设置窗口大小（像素）
tiledlayout(5, 3, 'Padding', 'compact', 'TileSpacing', 'compact');  % 5 行 3 列
mat_list = mat_list';  %

for i = 1:15
    nexttile(i);
    imagesc(mat_list{i}); % 重复使用上面的 mat_list
    xticks([]); yticks([]);
    caxis([0 1]);
    axis square;
end

colormap(othercolor('PuBu8'));
cb = colorbar;
cb.Layout.Tile = 'east';         % ← 这是关键，放在 tiledlayout 外侧
cb.FontSize = 10;
cb.FontWeight = 'bold';
caxis([0 1]);

exportgraphics(fig2, 'fig2.jpg', 'Resolution', 1500);