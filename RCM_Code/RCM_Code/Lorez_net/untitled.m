% p=squeeze(NCMindex(1,:,:))
% p=p(:)
figure
p=[squeeze(CCMindex(1,:,:));squeeze(CCMindex(2,:,:));squeeze(CCMindex(3,:,:));squeeze(CCMindex(4,:,:));squeeze(CCMindex(5,:,:))]
label=[connectivityMatrix{1};connectivityMatrix{2};connectivityMatrix{3};connectivityMatrix{4};connectivityMatrix{5}]
p=p(:)
label=label(:)
plot_roc(p,label,[0.1,0.7,0.9],'TE')
function  auc = plot_roc( predict, ground_truth ,color,name )
% INPUTS
%  predict       - 分类器对测试集的分类结果
%  ground_truth - 测试集的正确标签,这里只考虑二分类，即0和1
% OUTPUTS
%  auc            - 返回ROC曲线的曲线下的面积

%初始点为（1.0, 1.0）
x = 1.0;
y = 1.0;
%计算出ground_truth中正样本的数目pos_num和负样本的数目neg_num
pos_num = sum(ground_truth==1);
neg_num = sum(ground_truth==0);
%根据该数目可以计算出沿x轴或者y轴的步长
x_step = 1.0/neg_num;
y_step = 1.0/pos_num;
%首先对predict中的分类器输出值按照从小到大排列
[predict,index] = sort(predict);
ground_truth = ground_truth(index);
%对predict中的每个样本分别判断他们是FP或者是TP
%遍历ground_truth的元素，
%若ground_truth[i]=1,则TP减少了1，往y轴方向下降y_step
%若ground_truth[i]=0,则FP减少了1，往x轴方向下降x_step
X(1)=1;
Y(1)=1;
for i=1:length(ground_truth)
    if ground_truth(i) == 1
        y = y - y_step;
    else
        x = x - x_step;
    end
    X(i+1)=x;
    Y(i+1)=y;
end
%画出图像     
% plot(X,Y,'-ro','LineWidth',2,'MarkerSize',3);

plot(X, Y, 'LineWidth', 2, 'Color', color, 'LineStyle', '--');

% 设置 X 轴和 Y 轴范围及刻度
xlim([0 1]);
xticks(0:0.2:1);
ylim([0 1]);
yticks(0:0.2:1);

% 获取当前坐标轴对象
ax = gca;

% 设置只显示 X 轴和 Y 轴
ax.XAxis.Visible = 'on';       % 显示 X 轴
ax.YAxis.Visible = 'on';       % 显示 Y 轴
ax.Box = 'off';                % 隐藏顶部和右侧框线

% 设置字体属性
ax.FontWeight = 'bold';        % 字体加粗
ax.FontSize = 15;              % 字体大小
ax.FontName = 'Times New Roman'; % 设置字体为 Times New Roman

% 添加坐标轴标签
xlabel('False Positive Rate', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('True Positive Rate', 'FontSize', 15, 'FontWeight', 'bold');


%计算小矩形的面积,返回auc
auc = -trapz(X,Y);   
end