
function f = NCM3(data1,data2,dim,resSize,len,scale,Win_1,Win_2,Wres_1, norm1, norm2)

% 检查 norm1 和 norm2 的默认值
if nargin < 10  % 如果 norm1 未提供，设置为 true
    norm1 = true;
end
if nargin < 11  % 如果 norm2 未提供，设置为 true
    norm2 = true;
end

data=[data1;data2];
Data=data;

if norm1==true
    [Data, ps] = mapstd(data);% data normalization
end

% load Reservior_Parameter.mat
% resSize = 200; % number of RC neurons
% generate weight matrix
% sigma = 0.44;
% Win1 = -1 + 2*rand(resSize,dim);
% adj1 = zeros(resSize,dim);
% for m=1:resSize
%     for n=1:dim
%         if(rand(1,1)<sigma)  
%             adj1(m,n)=1;  
%         end
%     end
% end
% Win = adj1.*Win1;
% 
% arhow_r =0.79; % spectral radius
% d = 0.2; % sparsity
% k=round(d*resSize);
% adj2 = zeros(resSize,resSize);
% for i = 1:resSize
%     num = randperm(resSize,k);
%     for j = 1:k
%         adj2(i,num(j)) = 1;
%     end
% end
% Wres1 = -1 + 2*rand(resSize,resSize);
% Wres2 = adj2.*Wres1 ;
% SR = max(abs(eig(Wres2))) ;
% Wres = Wres2 .* ( arhow_r/SR);  

Win=Win_1;
Wres=Wres_1;
gamma = 0.44; % leaky rate
% gamma = 0.1; % leaky rate

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
dim=resSize;
resSize = 1000; % number of RC neurons

load("Embedding_reservior_parameter.mat")

Win = Win_2;
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

% vv = rtotal(:,len+1:len+testlen);
% testoutput = Wout*vv;
% 
% testdata = Data1(:,len+1:len+testlen);
% Mse2=mean(mean((testoutput-testdata).^2,2));
%     MyMse(runtime)=mean(sum((testoutput-testdata).^2,1))

% testdata_mean=mean(testdata,2);
% SUM=0;
% for i=1:testlen
%     SUM=SUM+(testdata(:,i)-testdata_mean).^2;
% end

% traindata_mean=mean(traindata,2);
% SUM=0;
% for i=1:trainlen
%     SUM=SUM+(traindata(:,i)-traindata_mean).^2;
% end
% 
% NMSE2=Mse1/mean(SUM);
% %     MyNMse(runtime)=MyMse(runtime)/sum(SUM);
% 
% INDEX=1/(1+100000*NMSE2);

dataset_mean = mean(traindata(:));
dataset_variance = mean((traindata(:) - dataset_mean).^2);
NMSE3=Mse1/ dataset_variance;
%     MyNMse(runtime)=MyMse(runtime)/sum(SUM);

INDEX=exp(-scale*NMSE3);

f=INDEX;
end

