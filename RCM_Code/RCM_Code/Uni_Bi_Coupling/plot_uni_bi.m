color_all = [120,149,193; 165,28,54]./255;
color_all = [1,0,0; 0,0,1];
i=2,j=1;
x=0:0.1:1

MEAN=(mean(INDEX,2))';
S = std(INDEX,0,2)'
yu=MEAN+S;
yl=MEAN-S;
color_mean = color_all(i,:);
color_between = color_all(i,:);
fig=figure('color',[1 1 1]);
plot(x,MEAN,'Linewidth',2,'DisplayName', 'RCM[X→Y]','color',color_mean);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between, 'linestyle', 'none', 'FaceAlpha',0.3); 
hold on;

MEAN=(mean(INDEX2,2))';
S = std(INDEX2,0,2)'
yu=MEAN+S;
yl=MEAN-S;
color_mean = color_all(i,:);
color_between = color_all(i,:);
plot(x,MEAN,'--','Linewidth',2,'DisplayName', 'RCM[Y→X]','color',color_mean);
axis([0,1,0,1])
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between, 'linestyle', 'none', 'FaceAlpha',0.3); 
hold on;

MEAN=(mean(x_y_MCM,2))';
S = std(x_y_RCM,0,2)'
yu=MEAN+S;
yl=MEAN-S;
color_mean = color_all(j,:);
color_between = color_all(j,:);
plot(x,MEAN,'Linewidth',2,'DisplayName', 'MCM[X→Y]','color',color_mean);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between, 'linestyle', 'none', 'FaceAlpha',0.3); 
hold on;

MEAN=(mean(y_x_MCM,2))';
S = std(y_x_RCM,0,2)'
yu=MEAN+S;
yl=MEAN-S;
color_mean = color_all(j,:);
color_between = color_all(j,:);
plot(x,MEAN,'--','Linewidth',2,'DisplayName', 'MCM[Y→X]','color',color_mean);
axis([0,1,0,1])
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between, 'linestyle', 'none', 'FaceAlpha',0.3); 

% 设置坐标轴属性
ax = gca;
ax.LineWidth = 1.3;
ax.Box = 'off';
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';
ax.FontName = 'Times New Roman';
ax.FontSize = 20;
ax.FontWeight = 'bold';      % 设置字体加粗

scatter(x, INDEX, 10, color_all(i,:),'filled');  % 将第 i 列的数据绘制成散点
scatter(x, INDEX2, 10, color_all(i,:),'filled');  % 将第 i 列的数据绘制成散点
scatter(x, x_y_MCM, 10, color_all(j,:),'filled');  % 将第 i 列的数据绘制成散点
scatter(x, y_x_MCM, 10, color_all(j,:),'filled');  % 将第 i 列的数据绘制成散点
gray_dash_line=[0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]
plot(x, gray_dash_line, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2)


% 添加图例
legend({'RCM[X→Y]','', 'RCM[Y→X]','', 'MCM[X→Y]','','MCM[Y→X]'}, ...
       'FontSize', 15, ...
       'Location', 'northwest', ...
       'TextColor', 'k', ...
       'Box', 'off', ...
       'Position',[0.124404761904762 0.725793657321778 0.296428565308451 0.23214285061473], ...
       'FontWeight', 'bold');  % 设置字体加粗
legend('show');  % 显示图例
exportgraphics(fig, 'fig2.jpg', 'Resolution', 1000);

%% 

load('rcm_no.mat')
MEAN=(mean(INDEX,2))';
S = std(INDEX,0,2)'
yu=MEAN+S;
yl=MEAN-S;
color_mean = color_all(i,:);
color_between = color_all(i,:);
figure('color',[1 1 1]);
plot(x,MEAN,'Linewidth',2,'DisplayName', 'CCM[X→Y]','color',color_mean);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between, 'linestyle', 'none', 'FaceAlpha',0.3); 
hold on;
INDEX_1=INDEX;

MEAN=(mean(INDEX2,2))';
S = std(INDEX2,0,2)'
yu=MEAN+S;
yl=MEAN-S;
color_mean = color_all(i,:);
color_between = color_all(i,:);
plot(x,MEAN,'--','Linewidth',2,'DisplayName', 'CCM[Y→X]','color',color_mean);
axis([0,0.45,0,1])
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between, 'linestyle', 'none', 'FaceAlpha',0.3); 
hold on;
INDEX_2=INDEX2;

load('mcm_no.mat')
MEAN=(mean(INDEX,2))';
S = std(INDEX,0,2)'
yu=MEAN+S;
yl=MEAN-S;
color_mean = color_all(j,:);
color_between = color_all(j,:);
plot(x,MEAN,'Linewidth',2,'DisplayName', 'CCM[X→Y]','color',color_mean);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between, 'linestyle', 'none', 'FaceAlpha',0.3); 
hold on;
INDEX_3=INDEX;

MEAN=(mean(INDEX2,2))';
S = std(INDEX2,0,2)'
yu=MEAN+S;
yl=MEAN-S;
color_mean = color_all(j,:);
color_between = color_all(j,:);
INDEX_4=INDEX2;
plot(x,MEAN,'--','Linewidth',2,'DisplayName', 'CCM[Y→X]','color',color_mean);
axis([0,0.45,0,1])
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between, 'linestyle', 'none', 'FaceAlpha',0.3); 

% 设置坐标轴属性
ax = gca;
ax.LineWidth = 1.3;
ax.Box = 'off';
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';
ax.FontName = 'Times New Roman';
ax.FontSize = 20;
ax.FontWeight = 'bold';      % 设置字体加粗

scatter(x, INDEX_1, 10, color_all(i,:),'filled');  % 将第 i 列的数据绘制成散点
scatter(x, INDEX_2, 10, color_all(i,:),'filled');  % 将第 i 列的数据绘制成散点
scatter(x, INDEX_3, 10, color_all(j,:),'filled');  % 将第 i 列的数据绘制成散点
scatter(x, INDEX_4, 10, color_all(j,:),'filled');  % 将第 i 列的数据绘制成散点

% 添加图例
legend({'RCM[X→Y]','', 'RCM[Y→X]','', 'MCM[X→Y]','','MCM[Y→X]'}, ...
       'FontSize', 15, ...
       'Location', 'northwest', ...
       'TextColor', 'k', ...
       'Box', 'off', ...
       'FontWeight', 'bold');  % 设置字体加粗
legend('show');  % 显示图例

% % 设置坐标轴属性
% ax = gca;
% ax.LineWidth = 1.3;
% ax.Box = 'off';
% ax.XAxis.Visible = 'on';
% ax.YAxis.Visible = 'on';
% ax.FontName = 'Times New Roman';
% ax.FontSize = 20;
% ax.FontWeight = 'bold';      % 设置字体加粗
