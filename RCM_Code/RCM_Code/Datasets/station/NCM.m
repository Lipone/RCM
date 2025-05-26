% with a standard non-delayed RC of dimension 200

function f = NCM(data1,data2,dim,len,scale, norm1, norm2)

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

load Reservior_Parameter.mat


gamma = 0.44; % leaky rate
% gamma = 0.1; % leaky rate

initialen = 1000;
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

initialen = 1000;%1000
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

