%robuteness plot
color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.905882352941177,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.305882352941177,0.541176470588235,0.317647058823529;
    1.0,0.5,0.0;
    0.5,0.7,0.4;
    0, 0, 0.75;];
CList = [153,153,253; 255,153,154]./255;
CList = [0, 0, 0.75; 1.0,0.5,0.0];

%% 画innner noise
load IN_NOISE_0_0.2.mat
figure('color',[1 1 1]);
% 第一组数据：NCM_values_for_ones
MEAN=zeros(1,11);
S=zeros(1,11);
for i = 1:11
    MEAN(:,i)=mean(all_ones{i});
    S(:,i)=std(all_ones{i});
end

x=0:0.02:0.2;
yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_ones = color_all(6,:);
color_between_ones = CList(2,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[X→Y]', 'color', color_mean_ones);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between_ones, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;

% 第二组数据：NCM_values_for_zeros
MEAN=zeros(1,11);
S=zeros(1,11);
for i = 1:11
    MEAN(:,i)=mean(all_zeros{i});
    S(:,i)=std(all_zeros{i});
end

yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_zeros = CList(1,:);
color_between_zeros = CList(1,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[Y→X]', 'color', color_mean_zeros);
fill([x fliplr(x)], [yu fliplr(yl)], color_between_zeros, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;
% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_zeros)
    y = [y, all_zeros{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, CList(1,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;

% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_ones)
    y = [y, all_ones{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, CList(2,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;
% 设置坐标轴属性
ax = gca;
ax.LineWidth = 1.3;
% xlabel('$\mathbf{\beta_{xy}}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
% 添加图例
legend({'RCM [with causation]','' , 'RCM [without causation]'}, 'Location', 'northeast', 'TextColor', 'k', 'Box', 'off','FontWeight', 'bold');  
legend('show');  % 隐藏顶部和右侧坐标轴
ax.Box = 'off';
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';
ax.FontName = 'Times New Roman';
ax.FontSize = 20;
ax.FontWeight = 'bold';      % 设置字体加粗
% ax.XTick=x
%% plot outter noise
load OUTER_NOISE_0-0.1-1(2).mat
figure('color',[1 1 1]);
% 第一组数据：NCM_values_for_ones
MEAN=zeros(1,11);
S=zeros(1,11);
for i = 1:11
    MEAN(:,i)=mean(all_ones{i});
    S(:,i)=std(all_ones{i});
end

x=0:0.1:1;
yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_ones = color_all(6,:);
color_between_ones = CList(2,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[X→Y]', 'color', color_mean_ones);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between_ones, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;

% 第二组数据：NCM_values_for_zeros
MEAN=zeros(1,11);
S=zeros(1,11);
for i = 1:11
    MEAN(:,i)=mean(all_zeros{i});
    S(:,i)=std(all_zeros{i});
end

yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_zeros = CList(1,:);
color_between_zeros = CList(1,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[Y→X]', 'color', color_mean_zeros);
fill([x fliplr(x)], [yu fliplr(yl)], color_between_zeros, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;
% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_zeros)
    y = [y, all_zeros{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, CList(1,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;

% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_ones)
    y = [y, all_ones{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, CList(2,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;
% 设置坐标轴属性
ax = gca;
ax.LineWidth = 1.3;
% xlabel('$\mathbf{\beta_{xy}}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
% 添加图例
legend({'RCM [with causation]','' , 'RCM [without causation]'}, 'Location', 'northeast', 'TextColor', 'k', 'Box', 'off','FontWeight', 'bold');  
legend('show'); 
% 隐藏顶部和右侧坐标轴
ax.Box = 'off';
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';
ax.FontName = 'Times New Roman';
ax.FontSize = 20;
ax.FontWeight = 'bold';      % 设置字体加粗

%% 画长度
load LEN_1000_10000.mat
figure('color',[1 1 1]);
% 第一组数据：NCM_values_for_ones
MEAN=zeros(1,10);
S=zeros(1,10);
for i = 1:10
    MEAN(:,i)=mean(all_ones{i});
    S(:,i)=std(all_ones{i});
end

x=1:1:10;
yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_ones = color_all(6,:);
color_between_ones = CList(2,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[X→Y]', 'color', color_mean_ones);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between_ones, 'linestyle', 'none', 'FaceAlpha', 0.3); 


% 第二组数据：NCM_values_for_zeros
MEAN=zeros(1,10);
S=zeros(1,10);
for i = 1:10
    MEAN(:,i)=mean(all_zeros{i});
    S(:,i)=std(all_zeros{i});
end

yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_zeros = CList(1,:);
color_between_zeros = CList(1,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[Y→X]', 'color', color_mean_zeros);
fill([x fliplr(x)], [yu fliplr(yl)], color_between_zeros, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;
hold on;
% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_zeros)
    y = [y, all_zeros{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, CList(1,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;

% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_ones)
    y = [y, all_ones{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, CList(2,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;

% 设置坐标轴属性
ax = gca;
ax.LineWidth = 1.3;
% xlabel('$\mathbf{\beta_{xy}}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
% 添加图例
legend({'RCM [with causation]','' , 'RCM [without causation]'}, 'Location', 'northeast', 'TextColor', 'k', 'Box', 'off','FontWeight', 'bold');  
legend('show');  
% 隐藏顶部和右侧坐标轴
ax.Box = 'off';
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';
ax.FontName = 'Times New Roman';
ax.FontSize = 20;
ax.XTick=x
ax.FontWeight = 'bold';      % 设置字体加粗

xlim([1 max(x)]);  % x 轴从 1 开始
%% plot emb
load EMB_10_20_50_100_200--1000.mat
all_ones=index_ones;
all_zeros=index_zeros;
figure('color',[1 1 1]);
% 第一组数据：NCM_values_for_ones
MEAN=zeros(1,10);
S=zeros(1,10);
for i = 1:10
    MEAN(:,i)=mean(all_ones{i});
    S(:,i)=std(all_ones{i});
end

x=1:1:10;
x_labels=[10,20,30,50,100,200,300,400,500,1000]
yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_ones = color_all(6,:);
color_between_ones = CList(2,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[X→Y]', 'color', color_mean_ones);
hold on;
fill([x fliplr(x)], [yu fliplr(yl)], color_between_ones, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;

% 第二组数据：NCM_values_for_zeros
MEAN=zeros(1,10);
S=zeros(1,10);
for i = 1:10
    MEAN(:,i)=mean(all_zeros{i});
    S(:,i)=std(all_zeros{i});
end

yu=MEAN+S;
yl=MEAN-S;
x_interp = linspace(min(x), max(x), 100);
yu_interp = interp1(x, yu, x_interp, 'spline');
yl_interp = interp1(x, yl, x_interp, 'spline');
color_mean_zeros = CList(1,:);
color_between_zeros = CList(1,:);
plot(x, MEAN, 'Linewidth', 2, 'DisplayName', 'CCM[Y→X]', 'color', color_mean_zeros);
fill([x fliplr(x)], [yu fliplr(yl)], color_between_zeros, 'linestyle', 'none', 'FaceAlpha', 0.3); 
hold on;
% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_zeros)
    y = [y, all_zeros{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, CList(1,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;

% 展开NCM_values_for_zeros中的数据，假设每个单元都是相同长度的数组
y = []; % 用来存储所有的 y 值

% 将NCM_values_for_zeros中每个数组的数据存储到 y 中
for i = 1:length(all_ones)
    y = [y, all_ones{i}]; % 将每个单元的内容添加到y数组
end
% 使用循环绘制每列数据的散点图
for i = 1:size(y', 2)  % 遍历每一列
    scatter(x, y(i, :), 10, CList(2,:),'filled');  % 将第 i 列的数据绘制成散点
end
hold on;
% 设置坐标轴属性
ax = gca;
ax.LineWidth = 1.3;
% xlabel('$\mathbf{\beta_{xy}}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
% 添加图例
legend({'RCM [with causation]','' , 'RCM [without causation]'}, 'Location', 'northeast', 'TextColor', 'k', 'Box', 'off','FontWeight', 'bold');  
legend('show');    
% 隐藏顶部和右侧坐标轴
ax.Box = 'off';
ax.XAxis.Visible = 'on';
ax.YAxis.Visible = 'on';
ax.FontName = 'Times New Roman';
ax.FontSize = 20;
ax.XTick=x
ax.XTickLabel = x_labels;
ax.XTickLabelRotation = 0; % 设置 X 轴刻度标签为水平显示
ax.FontWeight = 'bold';      % 设置字体加粗

xlim([1 max(x)]);  % x 轴从 1 开始