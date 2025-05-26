clear;
clc;
p=1;
k=5;
Thei=10;
numFiles = 11;  % 共11个文件 data1.mat 到 data11.mat
% 获取系统的物理核心数并启动并行池
numCores = feature('numcores');
parpool('local', numCores);  % 使用所有核心

for a = 1:numFiles
    filename = sprintf('Data_uni/data%d.mat', a);
    fileData = load(filename);
    dataCell = fileData.dataCell;
    parfor q=1:50
        y = dataCell{q};
        data=y(:,1000:6000);

        [x_y_RCM(a,q),y_x_RCM(a,q)]=RCM(data,length(data))

        x_y_MCM(a,q)=MCM(data(1,:), data(2,:), p, Thei)
        y_x_MCM(a,q)=MCM(data(2,:), data(1,:), p, Thei)

%         x_y_GC(a,q)=GCmy(data(1,:), data(2,:), p)
%         y_x_GC(a,q)=GCmy(data(2,:), data(1,:), p)
% 
%         x_y_TE(a,q)=TEmy(data(1,:), data(2,:), p,k,Thei)
%         y_x_TE(a,q)=TEmy(data(2,:), data(1,:), p,k,Thei)

    end

end

save('Uni2.mat', 'x_y_RCM', 'y_x_RCM', 'x_y_MCM', 'y_x_MCM');
%% 

for a = 1:numFiles
    filename = sprintf('Data_Bi/data%d.mat', a);
    fileData = load(filename);
    dataCell = fileData.dataCell;
    parfor q=1:50
        y = dataCell{q};
        data=y(:,1000:6000);

        [x_y_RCM(a,q),y_x_RCM(a,q)]=RCM(data,length(data))

%         x_y_MCM(a,q)=MCM(data(1,:), data(2,:), p, Thei)
%         y_x_MCM(a,q)=MCM(data(2,:), data(1,:), p, Thei)
% 
%         x_y_GC(a,q)=GCmy(data(1,:), data(2,:), p)
%         y_x_GC(a,q)=GCmy(data(2,:), data(1,:), p)
% 
%         x_y_TE(a,q)=TEmy(data(1,:), data(2,:), p,k,Thei)
%         y_x_TE(a,q)=TEmy(data(2,:), data(1,:), p,k,Thei)

    end

end

save('Bi3.mat', 'x_y_RCM', 'y_x_RCM');
