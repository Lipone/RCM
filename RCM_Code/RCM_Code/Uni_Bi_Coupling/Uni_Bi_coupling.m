%t实验好norm后进行画图
clear;
clc;


test_len=5000;
noise=0.015;
norm1=false;
norm2=true;

for a=1:11
    for q=1:50
        u12=0.1;
        u21=(a-1)*0.1;
        u31=0;
        u32=0;
        u13=0;
        u23=0;
        rx=3.8;
        ry=3.6;
        n=10000;
        y=rand(1,3)';

        for i=1:n
            true=1;

            y(1,i+1) = y(1,i)*(rx-(rx-u12-u13)*y(1,i)-u12*y(2,i)-u13*y(3,i))+normrnd(0,noise);
            y(2,i+1) = y(2,i)*(ry-(ry-u21-u23)*y(2,i)-u21*y(1,i)-u23*y(3,i))+normrnd(0,noise);

            if (y(1,i+1)<=100)&&(y(1,i+1)>-100)&&(y(2,i+1)<=100)&&(y(2,i+1)>-100)
                true=1;
            else true=0;
                while(true==0)
                    for i=1:n
                        y(1,i+1) = y(1,i)*(rx-(rx-u12-u13)*y(1,i)-u12*y(2,i)-u13*y(3,i))+normrnd(0,noise);
                        y(2,i+1) = y(2,i)*(ry-(ry-u21-u23)*y(2,i)-u21*y(1,i)-u23*y(3,i))+normrnd(0,noise);

                        if (y(1,i+1)<=100)&&(y(1,i+1)>-100)&&(y(2,i+1)<=100)&&(y(2,i+1)>-100)
                            true=1;
                        else true=0;break;
                        end
                    end
                end
                break;
            end
        end

        data=y(:,1000:10000);

        Data=data;
        if norm1==true
            [Data, ps] = mapstd(data);% data normalization
        end
        resSize = 200; % number of RC neurons
        res=200;
        dim = 1;        % 输入维度
        sigma = 0.44;    % 输入连接密度
        arhow_r = 0.67;  % 目标谱半径
        d = 0.05;        % 稀疏度

        %                 load Reservior_Parameter.mat

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

        Data1=rtrain1;
        Data2=rtrain2;
        if norm2==true
            [Data1, ps] = mapstd(rtrain1);% data normalization
            [Data2, ps] = mapstd(rtrain2);% data normalization
        end
        resSize = 1000; % number of RC neurons

        load("Embedding_reservior_parameter.mat")


        initialen = 100;%1000
        trainlen = test_len-initialen;%1000
        len = initialen+trainlen;
        testlen = 0000;%1000

        r = zeros(resSize,1);
        rtotal = zeros(resSize,len+testlen);
        gamma=1;

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

        beta =1e-5; % regularization parameter 1e-8
        Wout = ((rtrain*rtrain' + beta*eye(resSize)) \ (rtrain*traindata'))';
        trainoutput = Wout*rtrain;
        Mse1(a,q)=mean(mean((trainoutput-traindata).^2,2));
        dataset_mean = mean(traindata(:));
        dataset_variance = mean((traindata(:) - dataset_mean).^2);
        NMSE3(a,q)=Mse1(a, q) / dataset_variance;

        INDEX(a,q)=exp(-5*NMSE3(a,q))

        gamma=0.1;


        %         if norm2==true
        [Data1, ps] = mapstd(rtrain2);% data normalization
        [Data2, ps] = mapstd(rtrain1);% data normalization
        %         end



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
        if a>=6
            beta =  1e-1; % regularization parameter 1e-8
        else
            beta =  1e-5
        end
% beta =  1e-1
        Wout = ((rtrain*rtrain' + beta*eye(resSize)) \ (rtrain*traindata'))';
        trainoutput = Wout*rtrain;
        Mse1(a,q)=mean(mean((trainoutput-traindata).^2,2));
        dataset_mean = mean(traindata(:));
        dataset_variance = mean((traindata(:) - dataset_mean).^2);
        NMSE3(a,q)=Mse1(a, q) / dataset_variance;

        INDEX2(a,q)=exp(-5*NMSE3(a,q))

    end
end

figure;
plot(INDEX)
hold on;
plot(INDEX2)

%%

t = (0:0.05:0.45);


color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.905882352941177,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.305882352941177,0.541176470588235,0.317647058823529;
    1.0,0.5,0.0;
    0.5,0.7,0.4;];

% Define x-axis values
x = 0:0.05:0.45;

% Create figure
figure('color',[1 1 1]);

% Plot for INDEX
MEAN1 = mean(INDEX, 2)';
S1 = std(INDEX, 0, 2)';
yu1 = MEAN1 + S1;
yl1 = MEAN1 - S1;
color_mean1 = color_all(6, :);
color_between1 = color_all(6, :);

% Plot the mean line and the shaded area for INDEX
plot(x, MEAN1, 'Linewidth', 1.5, 'color', color_mean1);
hold on;
fill([x fliplr(x)], [yu1 fliplr(yl1)], color_between1, 'linestyle', 'none', 'FaceAlpha', 0.5);

% Plot for INDEX2
MEAN2 = mean(INDEX2, 2)';
S2 = std(INDEX2, 0, 2)';
yu2 = MEAN2 + S2;
yl2 = MEAN2 - S2;
color_mean2 = color_all(2, :);
color_between2 = color_all(2, :);

% Plot the mean line and the shaded area for INDEX2
plot(x, MEAN2, 'Linewidth', 1.5, 'color', color_mean2);
fill([x fliplr(x)], [yu2 fliplr(yl2)], color_between2, 'linestyle', 'none', 'FaceAlpha', 0.5);

% Set the axis limits
axis([0, 0.45, 0, 1]);

% Add legend
legend('INDEX', 'INDEX2');

% Add labels if needed
xlabel('X-axis label');
ylabel('Y-axis label');
title('Combined Plot for INDEX and INDEX2');
%%
clc;
clear;
test_len=5000;

n=5000;
rx=3.8;ry=3.6;rz=3.78;
beta_xy=0;beta_yz=0;beta_xz=0;beta_yx=0.1;
p=1;k=5;Thei=10;

for step=0:0.05:0.45
    beta_xy=step;
    for k=1:3
        true=1;
        x(k,1)=rand(1,1);y(k,1)=rand(1,1);z(k,1)=rand(1,1);

        for i=1:n

            x(k,i+1)=rx*x(k,i)*(1-(1-beta_yx/rx)*x(k,i)-beta_yx/rx*y(k,i))+normrnd(0,0.015);
            y(k,i+1)=ry*y(k,i)*(1-(1-beta_xy/ry)*y(k,i)-beta_xy/ry*x(k,i))+normrnd(0,0.015);
            % z(k,i+1)=rz*z(k,i)*(1-(1-(beta_xy+beta_yz)/rz)*z(k,i)-beta_xz/rz*x(k,i)-beta_yz/rz*y(k,i))+normrnd(0,0.01);

            if (x(k,i+1)<=100)&&(x(k,i+1)>-100)&&(y(k,i+1)<=100)&&(y(k,i+1)>-100)
                true=1;
            else true=0;
                while(true==0)
                    for i=1:n
                        % noise=randn(3,1)*0.02;
                        x(k,i+1)=rx*x(k,i)*(1-(1-beta_yx/rx)*x(k,i)-beta_yx/rx*y(k,i))+normrnd(0,0.015);
                        y(k,i+1)=ry*y(k,i)*(1-(1-beta_xy/ry)*y(k,i)-beta_xy/ry*x(k,i))+normrnd(0,0.015);
                        % z(k,i+1)=rz*z(k,i)*(1-(1-(beta_xy+beta_yz)/rz)*z(k,i)-beta_xz/rz*x(k,i)-beta_yz/rz*y(k,i))+normrnd(0,0.01);

                        if (x(k,i+1)<=100)&&(x(k,i+1)>-100)&&(y(k,i+1)<=100)&&(y(k,i+1)>-100)
                            true=1;
                        else true=0;break;
                        end
                    end
                end
                break;
            end

        end
        x_y_MCM(round(step/0.05+1),k)=MCM(x(k,1000:test_len), y(k,1000:test_len), p, Thei)
        y_x_MCM(round(step/0.05+1),k)=MCM(y(k,1000:test_len), x(k,1000:test_len), p, Thei)
        % x_y_EC(round(step/0.05+1),k)=ECmy(x(k,1000:5000), y(k,1000:5000), p, Thei)
        % y_x_EC(round(step/0.05+1),k)=ECmy(y(k,1000:5000), x(k,1000:5000), p, Thei)


    end
end

INDEX=x_y_MCM;
INDEX2=y_x_MCM;

t = (0:0.05:0.45);


color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.905882352941177,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.305882352941177,0.541176470588235,0.317647058823529;
    1.0,0.5,0.0;
    0.5,0.7,0.4;];

% Define x-axis values
x = 0:0.05:0.45;

% Create figure
figure('color',[1 1 1]);

% Plot for INDEX
MEAN1 = mean(INDEX, 2)';
S1 = std(INDEX, 0, 2)';
yu1 = MEAN1 + S1;
yl1 = MEAN1 - S1;
color_mean1 = color_all(6, :);
color_between1 = color_all(6, :);

% Plot the mean line and the shaded area for INDEX
plot(x, MEAN1, 'Linewidth', 1.5, 'color', color_mean1);
hold on;
fill([x fliplr(x)], [yu1 fliplr(yl1)], color_between1, 'linestyle', 'none', 'FaceAlpha', 0.5);

% Plot for INDEX2
MEAN2 = mean(INDEX2, 2)';
S2 = std(INDEX2, 0, 2)';
yu2 = MEAN2 + S2;
yl2 = MEAN2 - S2;
color_mean2 = color_all(2, :);
color_between2 = color_all(2, :);

% Plot the mean line and the shaded area for INDEX2
plot(x, MEAN2, 'Linewidth', 1.5, 'color', color_mean2);
fill([x fliplr(x)], [yu2 fliplr(yl2)], color_between2, 'linestyle', 'none', 'FaceAlpha', 0.5);

% Set the axis limits
axis([0, 0.45, 0, 1]);

% Add legend
legend('INDEX', 'INDEX2');

% Add labels if needed
xlabel('X-axis label');
ylabel('Y-axis label');
title('Combined Plot for INDEX and INDEX2');






