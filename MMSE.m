close all;
clear all;

rng('shuffle'); 
Nantennas = [5,10];
K = 4;
nbrOfMonteCarloRealizations = 2;
%Combined channel matrix will be (K x K*N). This matrix gives the
%normalized variance of each channel element
channelVariances = [1 1 1 1];
weights = [1 1 1 1]'; ones(K,1);

%Range of SNR values
PdB = -10:1:30; %dB scale
P = 10.^(PdB/10); %Linear scale
sumRateMMSE = zeros(length(P),nbrOfMonteCarloRealizations);
sumrateMMSE_G = zeros(length(P),nbrOfMonteCarloRealizations);
sumRateMMSE_10 = zeros(length(P),nbrOfMonteCarloRealizations);


%n = 5;
N = 5;
    
    %Pre-generation of Rayleigh fading channel realizations (unit variance)
    Hall = (randn(K,N,nbrOfMonteCarloRealizations)+1i*randn(K,N,nbrOfMonteCarloRealizations))/sqrt(2);
    
    %Go through all channel realizations
    for m = 1:nbrOfMonteCarloRealizations
        
    %Output the progress
    disp(['Progress: N = ' num2str(N) ', ' num2str(m) ' out of ' num2str(nbrOfMonteCarloRealizations) ' realizations.']);
        
        
        %Generate channel matrix for m:th realization
        H = repmat(sqrt(channelVariances)',[1 N]) .* Hall(:,:,m);
		
		
        %Go through all transmit powers
        for pind = 1:length(P)
            
            %Compute normalized beamforming vectors for transmit MMSE
            %beamforming (which is the same as regularized ZFBF and SLNR-MAX
            %beamforming). Note that it varies with the transmit power.
            wSLNRMAX = functionSLNRMAX(H,P(pind)*ones(K,1));
                 
            %Calculate power allocation with transmit MMSE beamforming (using Theorem 3.5 in [7])
            rhos = diag(abs(H*wSLNRMAX).^2)';
            powerAllocationwSLNRMAX_sumrate = functionHeuristicPowerAllocation(rhos,P(pind),weights);
            
            %Calculate sum rate with transmit MMSE beamforming
            W = kron(sqrt(powerAllocationwSLNRMAX_sumrate),ones(N,1)).*wSLNRMAX;
            channelGains = abs(H*W).^2;
            signalGains = diag(channelGains);
            interferenceGains = sum(channelGains,2)-signalGains;
            rates = log2(1+signalGains./(interferenceGains+1));
            sumRateMMSE(pind,m) = weights'*rates;
        end

		GH = load('channel_GA.mat');
		HG = GH.HG;
for pind = 1:length(P)
            
            %Compute normalized beamforming vectors for transmit MMSE
            %beamforming (which is the same as regularized ZFBF and SLNR-MAX
            %beamforming). Note that it varies with the transmit power.
            wSLNRMAX = functionSLNRMAX(HG,P(pind)*ones(K,1));
                 
            %Calculate power allocation with transmit MMSE beamforming (using Theorem 3.5 in [7])
            rhos = diag(abs(HG*wSLNRMAX).^2)';
            powerAllocationwSLNRMAX_sumrate = functionHeuristicPowerAllocation(rhos,P(pind),weights);
            
            %Calculate sum rate with transmit MMSE beamforming
            W = kron(sqrt(powerAllocationwSLNRMAX_sumrate),ones(N,1)).*wSLNRMAX;
            channelGains = abs(HG*W).^2;
            signalGains = diag(channelGains);
            interferenceGains = sum(channelGains,2)-signalGains;
            rates = log2(1+signalGains./(interferenceGains+1));
            sumrateMMSE_G(pind,m) = weights'*rates;
        end
		end
		
		
N = 10;
Hall = (randn(K,N,nbrOfMonteCarloRealizations)+1i*randn(K,N,nbrOfMonteCarloRealizations))/sqrt(2);
for m = 1:nbrOfMonteCarloRealizations
        
    %Output the progress
    disp(['Progress: N = ' num2str(N) ', ' num2str(m) ' out of ' num2str(nbrOfMonteCarloRealizations) ' realizations.']);
        
        
        %Generate channel matrix for m:th realization
        H = repmat(sqrt(channelVariances)',[1 N]) .* Hall(:,:,m);
		for pind = 1:length(P)
            
            %Compute normalized beamforming vectors for transmit MMSE
            %beamforming (which is the same as regularized ZFBF and SLNR-MAX
            %beamforming). Note that it varies with the transmit power.
            wSLNRMAX = functionSLNRMAX(H,P(pind)*ones(K,1));
                 
            %Calculate power allocation with transmit MMSE beamforming (using Theorem 3.5 in [7])
            rhos = diag(abs(H*wSLNRMAX).^2)';
            powerAllocationwSLNRMAX_sumrate = functionHeuristicPowerAllocation(rhos,P(pind),weights);
            
            %Calculate sum rate with transmit MMSE beamforming
            W = kron(sqrt(powerAllocationwSLNRMAX_sumrate),ones(N,1)).*wSLNRMAX;
            channelGains = abs(H*W).^2;
            signalGains = diag(channelGains);
            interferenceGains = sum(channelGains,2)-signalGains;
            rates = log2(1+signalGains./(interferenceGains+1));
            sumRateMMSE_10(pind,m) = weights'*rates;
        end
end
%%绘图部分


figure; hold on; box on;
plot(PdB,mean(sumRateMMSE(:,:),2),'r','LineWidth',1);
plot(PdB,mean(sumrateMMSE_G(:,:),2),'b--','LineWidth',1);
plot(PdB,mean(sumRateMMSE_10(:,:),2),'b','LineWidth',1);
legend('regular MMSE','MMSE with Antenna Selection','MMSE with 10 antennas','Location','NorthWest');
title(['A Comaprison of regular MMSE & MMSE with Antenna Selection']);

xlabel('Average SNR [dB]')
ylabel('Average Sum Rate [bit/channel use]');
