function S_rates = fun(xs,H)
%该函数用来计算个体的适应度值

%input ： x - 染色体
%input ： H - 实际信道参数

%output ：S_rates 不同个体对应的和速率

x = sort01(xs); %预处理
K = length(H(:,1));
num = length(H(1,:))/2;
H_selected = zeros(length(H(:,1)),num);
count = 0;
for i = 1:length(x)
	if x(i) > 0
	count = count+1;
	H_selected(:,count) = H(:,i);

	end
end

%%计算和速率
PdB = -10:1:30; %dB scale
P = 10.^(PdB/10); %Linear scale
weights = [1 1 1 1]'; ones(K,1);
for pind = 1:length(P)
wSLNRMAX = functionSLNRMAX(H_selected,P(pind)*ones(K,1));
%Calculate power allocation with transmit MMSE beamforming (using Theorem 3.5 in [7])
rhos = diag(abs(H_selected*wSLNRMAX).^2)';
powerAllocationwSLNRMAX_sumrate = functionHeuristicPowerAllocation(rhos,P(pind),weights);
            
%Calculate sum rate with transmit MMSE beamforming
W = kron(sqrt(powerAllocationwSLNRMAX_sumrate),ones(num,1)).*wSLNRMAX;
channelGains = abs(H_selected*W).^2;
signalGains = diag(channelGains);
interferenceGains = sum(channelGains,2)-signalGains;
rates = log2(1+signalGains./(interferenceGains+1));
S_rates = weights'*rates;
end