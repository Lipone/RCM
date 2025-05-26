u12=0.1;

rx=3.8;
ry=3.6;
n=10000;
noise=0.015;

% 结果存储
numGroups = 50;
dataCell = cell(numGroups, 1);
for a=1:11
    u21=(a-1)*0.1;
    for g = 1:numGroups
        valid = false;
        while ~valid
            y = zeros(2, n+1);
            y(:,1) = rand(2,1);  % 初始值
            valid = true;

            for i = 1:n
                y(1,i+1) = y(1,i)*(rx - (rx - u12 )*y(1,i) - u12*y(2,i)) + normrnd(0, noise);
                y(2,i+1) = y(2,i)*(ry - (ry - u21 )*y(2,i) - u21*y(1,i)) + normrnd(0, noise);

                if abs(y(1,i+1)) > 100 || abs(y(2,i+1)) > 100
                    valid = false;
                    break;
                end
            end
        end
        dataCell{g} = y;
    end
    % 构建保存路径并保存
    filename = sprintf('Data_Bi/data%d.mat', a);
    save(filename, 'dataCell');
    fprintf('已保存: %s (u21 = %.1f)\n', filename, u21);
end

