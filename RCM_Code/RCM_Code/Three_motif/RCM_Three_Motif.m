clear;
clc;

q1_range = 0.1:0.1:1;
q2_range = 0.1:0.1:1;

p=3;k=5;Thei=10;
couple_strength1=0.1
couple_strength2=0.1
trails=1
data_length=5000;

for q1 = 0.44
    for q2 = 0.9

        %% Fan-in
        couple_strength1=0.1
        couple_strength2=0.1
        beta=1e-8;q1 = 0.44;q2 = 0.9;
        for a=1:1
            for q=1:5
                u12=0;
                u21=couple_strength1;
                u31=0;
                u13=0;
                u23=couple_strength2;
                u32=0;
                rx=3.6;
                ry=3.72;
                rz=3.68;

                y=rand(1,3)';

                for i=1:data_length

                    y(1,i+1) = y(1,i)*(rx-(rx)*y(1,i)-u12*y(2,i)-u13*y(3,i))+normrnd(0,0.005);
                    y(2,i+1) = y(2,i)*(ry-(ry)*y(2,i)-u21*y(1,i)-u23*y(3,i))+normrnd(0,0.005);
                    y(3,i+1) = y(3,i)*(rz-(rz)*y(3,i)-u31*y(1,i)-u32*y(2,i))+normrnd(0,0.005);
                    %             y(1,i+1)=y(1,i)*(rx-rx*y(1,i)-u12*y(2,i));
                    %             y(2,i+1)=y(2,i)*(ry-ry*y(2,i)-u21*y(1,i));
                end

                data=y(:,1:data_length);
                %         CCM_XtoY(a,q)=MCM(data(1,1:5000),data(2,1:5000), p, Thei);
                %         CCM_YtoX(a,q)=MCM(data(2,1:5000),data(1,1:5000), p, Thei);
                %         CCM_XtoZ(a,q)=MCM(data(1,1:5000),data(3,1:5000), p, Thei);
                %         CCM_ZtoX(a,q)=MCM(data(3,1:5000),data(1,1:5000), p, Thei);
                %         CCM_YtoZ(a,q)=MCM(data(2,1:5000),data(3,1:5000), p, Thei);
                %         CCM_ZtoY(a,q)=MCM(data(3,1:5000),data(2,1:5000), p, Thei);

                RCM_XtoY(a,q)=RCM(data(1,1:data_length),data(2,1:data_length),data_length,q1,q2,beta);
                RCM_YtoX(a,q)=RCM(data(2,1:data_length),data(1,1:data_length),data_length,q1,q2,beta);
                RCM_XtoZ(a,q)=RCM(data(1,1:data_length),data(3,1:data_length),data_length,q1,q2,beta);
                RCM_ZtoX(a,q)=RCM(data(3,1:data_length),data(1,1:data_length),data_length,q1,q2,beta);
                RCM_YtoZ(a,q)=RCM(data(2,1:data_length),data(3,1:data_length),data_length,q1,q2,beta);
                RCM_ZtoY(a,q)=RCM(data(3,1:data_length),data(2,1:data_length),data_length,q1,q2,beta);
            end
        end

        % CCM_XtoY=mean(CCM_XtoY);
        % CCM_XtoZ=mean(CCM_XtoZ);
        % CCM_YtoX=mean(CCM_YtoX);
        % CCM_YtoZ=mean(CCM_YtoZ);
        % CCM_ZtoX=mean(CCM_ZtoX);
        % CCM_ZtoY=mean(CCM_ZtoY);

        RCM_XtoY=mean(RCM_XtoY);
        RCM_XtoZ=mean(RCM_XtoZ);
        RCM_YtoX=mean(RCM_YtoX);
        RCM_YtoZ=mean(RCM_YtoZ);
        RCM_ZtoX=mean(RCM_ZtoX);
        RCM_ZtoY=mean(RCM_ZtoY);

        RCM_Fan_in_matrix = [
            0,         RCM_XtoY, RCM_XtoZ;
            RCM_YtoX,  0,        RCM_YtoZ;
            RCM_ZtoX,  RCM_ZtoY, 0
            ];

        % CCM_Fan_in_matrix = [
        %     0,         CCM_XtoY, CCM_XtoZ;
        %     CCM_YtoX,  0,        CCM_YtoZ;
        %     CCM_ZtoX,  CCM_ZtoY, 0
        % ];
        m1=RCM_XtoY;
        m2=RCM_ZtoY;

        %% Fan-out

        for a=1:1
            for q=1:trails
                u12=0;
                u21=couple_strength1;
                u31=couple_strength2;
                u13=0;
                u23=0;
                u32=0;
                rx=3.6;
                ry=3.72;
                rz=3.68;

                y=rand(1,3)';

                for i=1:data_length

                    y(1,i+1) = y(1,i)*(rx-(rx)*y(1,i)-u12*y(2,i)-u13*y(3,i))+normrnd(0,0.005);
                    y(2,i+1) = y(2,i)*(ry-(ry)*y(2,i)-u21*y(1,i)-u23*y(3,i))+normrnd(0,0.005);
                    y(3,i+1) = y(3,i)*(rz-(rz)*y(3,i)-u31*y(1,i)-u32*y(2,i))+normrnd(0,0.005);
                    %             y(1,i+1)=y(1,i)*(rx-rx*y(1,i)-u12*y(2,i));
                    %             y(2,i+1)=y(2,i)*(ry-ry*y(2,i)-u21*y(1,i));
                end

                data=y(:,1:data_length);
                %         CCM_XtoY(a,q)=MCM(data(1,1:5000),data(2,1:5000), p, Thei);
                %         CCM_YtoX(a,q)=MCM(data(2,1:5000),data(1,1:5000), p, Thei);
                %         CCM_XtoZ(a,q)=MCM(data(1,1:5000),data(3,1:5000), p, Thei);
                %         CCM_ZtoX(a,q)=MCM(data(3,1:5000),data(1,1:5000), p, Thei);
                %         CCM_YtoZ(a,q)=MCM(data(2,1:5000),data(3,1:5000), p, Thei);
                %         CCM_ZtoY(a,q)=MCM(data(3,1:5000),data(2,1:5000), p, Thei);

                RCM_XtoY(a,q)=RCM(data(1,1:data_length),data(2,1:data_length),data_length,q1,q2,beta);
                RCM_YtoX(a,q)=RCM(data(2,1:data_length),data(1,1:data_length),data_length,q1,q2,beta);
                RCM_XtoZ(a,q)=RCM(data(1,1:data_length),data(3,1:data_length),data_length,q1,q2,beta);
                RCM_ZtoX(a,q)=RCM(data(3,1:data_length),data(1,1:data_length),data_length,q1,q2,beta);
                RCM_YtoZ(a,q)=RCM(data(2,1:data_length),data(3,1:data_length),data_length,q1,q2,beta);
                RCM_ZtoY(a,q)=RCM(data(3,1:data_length),data(2,1:data_length),data_length,q1,q2,beta);
            end
        end

        % CCM_XtoY=mean(CCM_XtoY);
        % CCM_XtoZ=mean(CCM_XtoZ);
        % CCM_YtoX=mean(CCM_YtoX);
        % CCM_YtoZ=mean(CCM_YtoZ);
        % CCM_ZtoX=mean(CCM_ZtoX);
        % CCM_ZtoY=mean(CCM_ZtoY);

        RCM_XtoY=mean(RCM_XtoY);
        RCM_XtoZ=mean(RCM_XtoZ);
        RCM_YtoX=mean(RCM_YtoX);
        RCM_YtoZ=mean(RCM_YtoZ);
        RCM_ZtoX=mean(RCM_ZtoX);
        RCM_ZtoY=mean(RCM_ZtoY);

        RCM_Fan_out_matrix = [
            0,         RCM_XtoY, RCM_XtoZ;
            RCM_YtoX,  0,        RCM_YtoZ;
            RCM_ZtoX,  RCM_ZtoY, 0
            ];

        % CCM_Fan_out_matrix = [
        %     0,         CCM_XtoY, CCM_XtoZ;
        %     CCM_YtoX,  0,        CCM_YtoZ;
        %     CCM_ZtoX,  CCM_ZtoY, 0
        % ];
        m3=RCM_XtoY;
        m4=RCM_XtoZ;

        %% Cascade
        beta=1e-8;
        for a=1:1
            for q=1:5
                u12=0;
                u21=0;
                u31=couple_strength1;
                u13=0;
                u23=couple_strength2;
                u32=0;
                rx=3.6;
                ry=3.72;
                rz=3.68;

                y=rand(1,3)';

                for i=1:data_length

                    y(1,i+1) = y(1,i)*(rx-(rx)*y(1,i)-u12*y(2,i)-u13*y(3,i))+normrnd(0,0.005);
                    y(2,i+1) = y(2,i)*(ry-(ry)*y(2,i)-u21*y(1,i)-u23*y(3,i))+normrnd(0,0.005);
                    y(3,i+1) = y(3,i)*(rz-(rz)*y(3,i)-u31*y(1,i)-u32*y(2,i))+normrnd(0,0.005);
                    %             y(1,i+1)=y(1,i)*(rx-rx*y(1,i)-u12*y(2,i));
                    %             y(2,i+1)=y(2,i)*(ry-ry*y(2,i)-u21*y(1,i));
                end

                data=y(:,1:data_length);
                %         CCM_XtoY(a,q)=MCM(data(1,1:5000),data(2,1:5000), p, Thei);
                %         CCM_YtoX(a,q)=MCM(data(2,1:5000),data(1,1:5000), p, Thei);
                %         CCM_XtoZ(a,q)=MCM(data(1,1:5000),data(3,1:5000), p, Thei);
                %         CCM_ZtoX(a,q)=MCM(data(3,1:5000),data(1,1:5000), p, Thei);
                %         CCM_YtoZ(a,q)=MCM(data(2,1:5000),data(3,1:5000), p, Thei);
                %         CCM_ZtoY(a,q)=MCM(data(3,1:5000),data(2,1:5000), p, Thei);

                RCM_XtoY(a,q)=RCM(data(1,1:data_length),data(2,1:data_length),data_length,q1,q2,beta);
                RCM_YtoX(a,q)=RCM(data(2,1:data_length),data(1,1:data_length),data_length,q1,q2,beta);
                RCM_XtoZ(a,q)=RCM(data(1,1:data_length),data(3,1:data_length),data_length,q1,q2,beta);
                RCM_ZtoX(a,q)=RCM(data(3,1:data_length),data(1,1:data_length),data_length,q1,q2,beta);
                RCM_YtoZ(a,q)=RCM(data(2,1:data_length),data(3,1:data_length),data_length,q1,q2,beta);
                RCM_ZtoY(a,q)=RCM(data(3,1:data_length),data(2,1:data_length),data_length,q1,q2,beta);
            end
        end

        % CCM_XtoY=mean(CCM_XtoY);
        % CCM_XtoZ=mean(CCM_XtoZ);
        % CCM_YtoX=mean(CCM_YtoX);
        % CCM_YtoZ=mean(CCM_YtoZ);
        % CCM_ZtoX=mean(CCM_ZtoX);
        % CCM_ZtoY=mean(CCM_ZtoY);

        RCM_XtoY=mean(RCM_XtoY);
        RCM_XtoZ=mean(RCM_XtoZ);
        RCM_YtoX=mean(RCM_YtoX);
        RCM_YtoZ=mean(RCM_YtoZ);
        RCM_ZtoX=mean(RCM_ZtoX);
        RCM_ZtoY=mean(RCM_ZtoY);

        RCM_Cascade_matrix = [
            0,         RCM_XtoY, RCM_XtoZ;
            RCM_YtoX,  0,        RCM_YtoZ;
            RCM_ZtoX,  RCM_ZtoY, 0
            ];

        % CCM_Cascade_matrix = [
        %     0,         CCM_XtoY, CCM_XtoZ;
        %     CCM_YtoX,  0,        CCM_YtoZ;
        %     CCM_ZtoX,  CCM_ZtoY, 0
        % ];
        %
        m5=RCM_XtoZ;
        m6=RCM_ZtoY;

        Cuasal_matrix=[RCM_Cascade_matrix;RCM_Fan_in_matrix;RCM_Fan_out_matrix];
        values = Cuasal_matrix(:);
        values = values(values > 0);  % 排除对角线的0
        top_six = maxk(values,6);
        if (ismember(m1, top_six) && ismember(m2, top_six) && ismember(m3, top_six) && ismember(m4, top_six) && ismember(m5, top_six) && ismember(m6, top_six))
            fprintf('q1 = %.2f, q2 = %.2f\n', q1, q2);
            disp('对应的RCM因果矩阵为：');
            disp(Cuasal_matrix);
        end

    end
end


