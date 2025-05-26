% 参数定义
resSize=[10,20,30,50,100,200,300,400,500,1000];      % reservoir size（自行设定）
dim = 3;        % 输入维度


% 生成并保存10组
for res=resSize
    for idx = 1:10
        sigma = 0.44;    % 输入连接密度
        arhow_r = 0.67;  % 目标谱半径
        d = 0.05;        % 稀疏度
        %% === 生成 Win ===
        Win1 = -1 + 2*rand(res,dim);
        adj1 = zeros(res,dim);
        for m=1:res
            for n=1:dim
                if(rand(1,1)<sigma)
                    adj1(m,n)=1;
                end
            end
        end
        Win = adj1.*Win1;                 % 应用掩码

        %% === 生成 Wres（内部连接）===
        k = round(d * res);                     % 每行保留 k 个非零
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

        dim2=res;

        sigma=0.01;
        Win1 = -1 + 2*rand(1000,dim2);
        adj1 = zeros(1000,dim2);
        for m=1:1000
            for n=1:dim2
                if(rand(1,1)<sigma)
                    adj1(m,n)=1;
                end
            end
        end
        Win_predict = adj1.*Win1;

        %% === 保存 ===
        folder = 'Weights';  % 文件夹名称
        if ~exist(folder, 'dir')
            mkdir(folder);   % 如果不存在，则创建该文件夹
        end

        filename = sprintf('ESN_weights_EmbSize_%d_0%d.mat', res,idx);
        save(fullfile(folder, filename), 'Win', 'Wres', 'Win_predict');

        fprintf('Saved: %s\n', filename);
    end
end

