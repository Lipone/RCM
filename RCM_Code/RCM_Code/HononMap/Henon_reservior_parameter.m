% reservoir parameter
arhow_r =0.67; % spectral radius
d = 0.05; % sparsity
resSize=200;
k = round(d*resSize);
inSize = 2; outSize = 2;
gamma = 0.1; % leaky rate
sigma = 0.44;

% generate weight matrix
Win1 = -1 + 2*rand(resSize,inSize);
adj1 = zeros(resSize,inSize);
for m=1:resSize
    for n=1:inSize
        if(rand(1,1)<sigma)  
            adj1(m,n)=1;  
        end
    end
end
Win = adj1.*Win1;


adj2 = zeros(resSize,resSize);

for i = 1:resSize
    num = randperm(resSize,k);
    for j = 1:k
        adj2(i,num(j)) = 1;
    end
end
Wres1 = -1 + 2*rand(resSize,resSize); 
Wres2 = adj2.*Wres1 ;
SR = max(abs(eig(Wres2))) ;
Wres = Wres2 .* ( arhow_r/SR);  