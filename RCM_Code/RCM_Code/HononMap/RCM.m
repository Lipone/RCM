% with a standard non-delayed RC of dimension 200
function f = RCM(data1,data2,dim,len,scale, norm1, norm2, Win, Wres)

% 检查 norm1 和 norm2 的默认值
if nargin < 6  % 如果 norm1 未提供，设置为 true
    norm1 = true;
end
if nargin < 7  % 如果 norm2 未提供，设置为 true
    norm2 = true;
end

data=[data1;data2];
Data=data;

if norm1==true
    [Data, ps] = mapstd(data);% data normalization
end

resSize = 200; % number of RC neurons

% if dim==1
%     load Reservior_Parameter.mat
% else
%     load NondelayedParameter.mat
% end
res=resSize;
sigma = 0.44;
Win=Win;
Wres=Wres;
% Win1 = -1 + 2*rand(res,dim);
% adj1 = zeros(res,dim);
% for m=1:res
%     for n=1:dim
%         if(rand(1,1)<sigma)  
%             adj1(m,n)=1;  
%         end
%     end
% end
% Win = adj1.*Win1;
% 
% %         arhow_r =0.79; % spectral radius
% arhow_r=0.67;
% d = 0.05; % sparsity
% k=round(d*res);
% adj2 = zeros(res,res);
% for i = 1:res
%     num = randperm(res,k);
%     for j = 1:k
%         adj2(i,num(j)) = 1;
%     end
% end
% Wres1 = -1 + 2*rand(res,res);
% Wres2 = adj2.*Wres1 ;
% SR = max(abs(eig(Wres2))) ;
% Wres = Wres2 .* ( arhow_r/SR);

gamma = 0.44; % leaky rate


initialen = 100;
trainlen = len-initialen;
len = initialen+trainlen;
testlen = 0;

r = zeros(resSize,1);
rtotal1 = zeros(resSize,len+testlen);
rtotal2 = zeros(resSize,len+testlen);

% training period
for i = 1:len+testlen
    ut = Data(1:dim,i);
    r = (1-gamma)*r + gamma*(tanh( Win*ut + Wres*r));
    rtotal1(:,i) = r;
end

r = zeros(resSize,1);
for i = 1:len+testlen
    ut = Data(dim+1:2*dim,i);
    r = (1-gamma)*r + gamma*(tanh( Win*ut + Wres*r));
    rtotal2(:,i) = r;
end


rtotal1(2:2:end,:)=rtotal1(2:2:end,:).^2; % half neurons are nonlinear(even terms)
rtotal2(2:2:end,:)=rtotal2(2:2:end,:).^2; % half neurons are nonlinear(even terms)
rtrain1 = rtotal1(:,initialen+1:len);
rtrain2 = rtotal2(:,initialen+1:len);

Data1=rtotal1;
Data2=rtotal2;

if norm2==true
    [Data1, ps] = mapstd(rtotal1);% data normalization
    [Data2, ps] = mapstd(rtotal2);% data normalization
end

% Data1=rtrain1;
% Data2=rtrain2;

resSize = 1000; % number of RC neurons

load("Embedding_reservior_parameter.mat")
gamma = 0.9;

initialen = 100;%1000
trainlen =len-initialen;%1000
len = initialen+trainlen;
testlen = 0;%1000

r = zeros(resSize,1);
rtotal = zeros(resSize,len+testlen);

% training period
for i = 1:len+testlen
    ut = Data2(:,i);
    r = (1-gamma)*r + gamma*(tanh( Win*ut + Wres*r));
    rtotal(:,i) = r;
end

rtotal(2:2:end,:)=rtotal(2:2:end,:).^2; % half neurons are nonlinear(even terms)

rtrain = rtotal(:,initialen+1:len);

% Tikhonov regularization to solve Wout
traindata =Data1(:,initialen+1:len) ;

beta = 1e-5; % regularization parameter 1e-8
Wout = ((rtrain*rtrain' + beta*eye(resSize)) \ (rtrain*traindata'))';
trainoutput = Wout*rtrain;
Mse1=mean(mean((trainoutput-traindata).^2,2));
dataset_mean = mean(traindata(:));
dataset_variance = mean((traindata(:) - dataset_mean).^2);
NMSE3=Mse1/ dataset_variance;
%     MyNMse(runtime)=MyMse(runtime)/sum(SUM);

INDEX=exp(-scale*NMSE3);

% 定义真实值矩阵和预测值矩阵，大小为 num x dim
y_true = traindata';
y_pred = trainoutput';

% 样本数和维度数
[N, D] = size(y_true);

% 初始化一个向量来存储每个样本的 NMSE
nmse_values = zeros(N, 1);

for i = 1:N
    % 获取第 i 行的真实值和预测值
    y_true_i = y_true(i, :);
    y_pred_i = y_pred(i, :);
    
    % 计算均方误差（MSE）
    mse = mean((y_true_i - y_pred_i).^2);
    
    % 计算真实值的方差
    variance = var(y_true_i, 1);
    
    % 计算归一化均方误差（NMSE）
    nmse_values(i) = mse / variance;
end
% mean(nmse_values)

f=exp(-scale*mean(nmse_values));
f=INDEX;
end




