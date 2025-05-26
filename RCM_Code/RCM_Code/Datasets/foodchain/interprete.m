%interprete
clear
clc
load ("foodchain.mat")
Data= data_FC(:,:)';
for i=1:4
    x = 1:794;
    y = Data(i,:);
    x1 = 1:0.1:794;    %x每改变0.1要插一个值
    % plot(x,y)
    % hold on
    % subplot(2,2,1)
    data(i,:) = interp1(x,y,x1,'spline');    %插值函数interp1:计算出这些插值点在y方向上的值存入y1中,y1很长
    % plot(x,y,x1,y1)   %完成插值后的曲线
    % xlabel('interp1(x,y,x1,''linear'')','color','b')
    % legend('y','y1','location','northwest')
    % title('method = ''linear''','color','r','fontsize',14)
end
save("processedData.mat","data")
% subplot(2,2,2)
% y1 = interp1(x,y,x1,'nearest');    %插值函数interp1:计算出这些插值点在y方向上的值存入y1中,y1很长
% plot(x,y,x1,y1)   %完成插值后的曲线
% xlabel('interp1(x,y,x1,''nearest'')','color','b')
% legend('y','y1','location','northwest')
% title('method = ''nearest''','color','r','fontsize',14)
%
% subplot(2,2,3)
% y1 = interp1(x,y,x1,'pchip');    %插值函数interp1:计算出这些插值点在y方向上的值存入y1中,y1很长
% plot(x,y,x1,y1)   %完成插值后的曲线
% xlabel('interp1(x,y,x1,''pchip'')','color','b')
% legend('y','y1','location','northwest')
% title('method = ''pchip''','color','r','fontsize',14)
%
% subplot(2,2,4)
% y1 = interp1(x,y,x1,'spline');    %插值函数interp1:计算出这些插值点在y方向上的值存入y1中,y1很长
% plot(x,y,x1,y1)   %完成插值后的曲线
% xlabel('interp1(x,y,x1,''spline'')','color','b')
% legend('y','y1','location','northwest')
% title('method = ''spline''','color','r','fontsize',14)
