% 参数定义
res = 200;       % reservoir size（自行设定）
dim = 3;        % 输入维度
sigma = 0.44;    % 输入连接密度
arhow_r = 0.67;  % 目标谱半径
d = 0.05;        % 稀疏度

% 生成并保存10组
for idx = 1:10
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

    %% === 保存 ===
    folder = 'Weights';  % 文件夹名称
    if ~exist(folder, 'dir')
        mkdir(folder);   % 如果不存在，则创建该文件夹
    end

    filename = sprintf('ESN_weights_%02d.mat', idx);
    save(fullfile(folder, filename), 'Win', 'Wres');

    fprintf('Saved: %s\n', filename);
end
