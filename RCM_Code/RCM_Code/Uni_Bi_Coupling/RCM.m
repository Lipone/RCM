% with a standard non-delayed RC of dimension 200
function [INDEX,INDEX2] = RCM(data,test_len)


Data=data;

resSize = 200; % number of RC neurons
res=200;
dim = 1;        % 输入维度
sigma = 0.44;    % 输入连接密度
arhow_r = 0.67;  % 目标谱半径
d = 0.05;        % 稀疏度

Win1 = -1 + 2*rand(res,dim);
adj1 = zeros(res,dim);
for m=1:res
    for n=1:dim
        if(rand(1,1)<sigma)
            adj1(m,n)=1;
        end
    end
end
Win = adj1.*Win1;
k = round(d * res);
adj2 = zeros(res,res);
for i = 1:res
    num = randperm(res,k);
    for j = 1:k
        adj2(i,num(j)) = 1;
    end
end
Wres1 = -1 + 2*rand(res,res);
Wres2 = adj2.*Wres1 ;
SR = max(abs(eig(Wres2))) ;
Wres = Wres2 .* ( arhow_r/SR);

initialen = 0000;
trainlen =  test_len;%6000
len = initialen+trainlen;
testlen = 0000;

r = zeros(resSize,1);
rtotal = zeros(resSize,len+testlen);
rtotal1 = zeros(resSize,len+testlen);
rtotal2 = zeros(resSize,len+testlen);

gamma=1;

% training period
for i = 1:len+testlen
    ut = Data(1,i);
    r = (1-gamma)*r + gamma*(tanh( Win*ut + Wres*r));
    rtotal1(:,i) = r;
end

r = zeros(resSize,1);
for i = 1:len+testlen
    ut = Data(2,i);
    r = (1-gamma)*r + gamma*(tanh( Win*ut + Wres*r));
    rtotal2(:,i) = r;
end


rtotal1(2:2:end,:)=rtotal1(2:2:end,:).^2; % half neurons are nonlinear(even terms)
rtotal2(2:2:end,:)=rtotal2(2:2:end,:).^2; % half neurons are nonlinear(even terms)
rtrain1 = rtotal1(:,initialen+1:len);
rtrain2 = rtotal2(:,initialen+1:len);


% Tikhonov regularization to solve Wout
% traindata =rtrain1 ;
% rtrain=rtrain2;
% beta = 1e-5; % regularization parameter
gamma=0.9;
Data1=rtrain1;
Data2=rtrain2;

[Data1, ps] = mapstd(rtrain1);% data normalization
[Data2, ps] = mapstd(rtrain2);% data normalization

resSize = 1000; % number of RC neurons

load("Embedding_reservior_parameter.mat")


initialen = 100;%1000
trainlen = test_len-initialen;%1000
len = initialen+trainlen;
testlen = 0000;%1000

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
NMSE3=Mse1 / dataset_variance;

INDEX=exp(-5*NMSE3)





[Data1, ps] = mapstd(rtrain2);% data normalization
[Data2, ps] = mapstd(rtrain1);% data normalization



r = zeros(resSize,1);
rtotal = zeros(resSize,len+testlen);
gamma=0.1;
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

beta =1; % regularization parameter 1e-8
Wout = ((rtrain*rtrain' + beta*eye(resSize)) \ (rtrain*traindata'))';
trainoutput = Wout*rtrain;
Mse1=mean(mean((trainoutput-traindata).^2,2));
dataset_mean = mean(traindata(:));
dataset_variance = mean((traindata(:) - dataset_mean).^2);
NMSE3=Mse1 / dataset_variance;

INDEX2=exp(-5*NMSE3)

end




