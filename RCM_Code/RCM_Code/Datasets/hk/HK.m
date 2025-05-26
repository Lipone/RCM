
clear
load ("HKdata.mat")

data=HKdata';
Data= data;
% [Data, ps] = mapstd(data);% data normalization


for L=0:100:1000
    initialen = L;
    trainlen =  1032-L;%6000 5156
    len = initialen+trainlen;
    testlen = 0;

    

    for m=1:5
        for n=1:5
            % training period
            resSize = 200; % number of RC neurons
            r = zeros(resSize,1);
    rtotal = zeros(resSize,len+testlen);
    rtotal1 = zeros(resSize,len+testlen);
    rtotal2 = zeros(resSize,len+testlen);
    load Reservior_Parameter.mat
            for i = 1:len+testlen
                ut = Data(m,i);
                r = (1-gamma)*r + gamma*(tanh( Win*ut + Wres*r));
                rtotal1(:,i) = r;
            end

            r = zeros(resSize,1);
            for i = 1:len+testlen
                ut = Data(n,i);
                r = (1-gamma)*r + gamma*(tanh( Win*ut + Wres*r));
                rtotal2(:,i) = r;
            end


            rtotal1(2:2:end,:)=rtotal1(2:2:end,:).^2; % half neurons are nonlinear(even terms)
            rtotal2(2:2:end,:)=rtotal2(2:2:end,:).^2; % half neurons are nonlinear(even terms)
            rtrain1 = rtotal1(:,initialen+1:len);
            rtrain2 = rtotal2(:,initialen+1:len);


            % Tikhonov regularization to solve Wout
            traindata =rtrain1 ;
            rtrain=rtrain2;
            beta = 1e-5; % regularization parameter
            Wout = ((rtrain*rtrain' + beta*eye(resSize)) \ (rtrain*traindata'))';
            trainoutput = Wout*rtrain2;
            mse1=mean(mean((trainoutput-traindata).^2,2));

            vv = rtotal2(:,len+1:len+testlen);
            testoutput = Wout*vv;

            testdata = rtotal1(:,len+1:len+testlen);
            mse2=mean(mean((testoutput-testdata).^2,2));

            Data1=rtotal1;
            Data2=rtotal2;
% 
%             [Data1, ps] = mapstd(rtotal1);% data normalization
%             [Data2, ps] = mapstd(rtotal2);% data normalization
            resSize = 1000; % number of RC neurons

            load("Embedding_reservior_parameter.mat")


            initialen = L;%1000
            trainlen = 1032-L;%1000
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
            Mse1(m,n,(L+100)/100)=mean(mean((trainoutput-traindata).^2,2));



            traindata_mean=mean(traindata,2);
            SUM=0;
            for i=1:trainlen
                SUM=SUM+(traindata(:,i)-traindata_mean).^2;
            end
            NMSE2=Mse1(m,n,(L+100)/100)/mean(SUM);
            %     MyNMse(runtime)=MyMse(runtime)/sum(SUM);

            INDEX(m,n,(L+100)/100)=1/(1+50000*NMSE2);
        end
    end

end

