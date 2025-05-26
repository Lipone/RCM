function out = RCM(data1,data2,data_length,q1,q2,beta)
Data = [data1; data2];  
x_original = 1:data_length;
x_interp = linspace(1, data_length, data_length * 2);  

Data_interp = zeros(2, length(x_interp));  

for i = 1:2
    Data_interp(i, :) = interp1(x_original, Data(i, :), x_interp, 'linear');
end

Data=Data_interp;


resSize = 200; % number of RC neurons

load Reservior_Parameter.mat

initialen = 100;
trainlen =  data_length-100;%6000
len = initialen+trainlen;
testlen = 0;

r = zeros(resSize,1);
rtotal = zeros(resSize,len+testlen);
rtotal1 = zeros(resSize,len+testlen);
rtotal2 = zeros(resSize,len+testlen);
g=q1;

% training period
for i = 1:len+testlen
    ut = Data(1,i);
    r = (1-g)*r + g*(tanh( Win*ut + Wres*r));
    rtotal1(:,i) = r;
end

r = zeros(resSize,1);
for i = 1:len+testlen
    ut = Data(2,i);
    r = (1-g)*r + g*(tanh( Win*ut + Wres*r));
    rtotal2(:,i) = r;
end


rtotal1(2:2:end,:)=rtotal1(2:2:end,:).^2; % half neurons are nonlinear(even terms)
rtotal2(2:2:end,:)=rtotal2(2:2:end,:).^2; % half neurons are nonlinear(even terms)
rtrain1 = rtotal1(:,initialen+1:len);
rtrain2 = rtotal2(:,initialen+1:len);


% Tikhonov regularization to solve Wout
traindata =rtrain1 ;
rtrain=rtrain2;
beta =beta ; % regularization parameter
Wout = ((rtrain*rtrain' + beta*eye(resSize)) \ (rtrain*traindata'))';
trainoutput = Wout*rtrain2;


[Data1, ps] = mapstd(rtrain1);% data normalization
[Data2, ps] = mapstd(rtrain2);% data normalization
resSize = 1000; % number of RC neurons

load("Embedding_reservior_parameter.mat")
g=q2;


initialen = 100;%1000
trainlen = data_length-200;%1000
len = initialen+trainlen;
testlen = 0;%1000

r = zeros(resSize,1);
rtotal = zeros(resSize,len+testlen);

% training period
for i = 1:len+testlen
    ut = Data2(:,i);
    r = (1-g)*r + g*(tanh( Win*ut + Wres*r));
    rtotal(:,i) = r;
end

rtotal(2:2:end,:)=rtotal(2:2:end,:).^2; % half neurons are nonlinear(even terms)

rtrain = rtotal(:,initialen+1:len);

% Tikhonov regularization to solve Wout
traindata =Data1(:,initialen+1:len) ;

beta =beta; % regularization parameter 1e-8
Wout = ((rtrain*rtrain' + beta*eye(resSize)) \ (rtrain*traindata'))';
trainoutput = Wout*rtrain;
Mse1=mean(mean((trainoutput-traindata).^2,2));
traindata_mean=mean(traindata,2);
dataset_mean = mean(traindata(:));
dataset_variance = mean((traindata(:) - dataset_mean).^2);
NMSE3=Mse1 / dataset_variance;

INDEX_new=exp(-5*NMSE3);

out=INDEX_new;
end