close all;
clear all;
%%
%暂时可用的版本1【
%
%
%%%
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

sumrateMMSE_G = zeros(length(P),nbrOfMonteCarloRealizations);
% N = 10;
N = 10;
    Hall = (randn(K,N,nbrOfMonteCarloRealizations)+1i*randn(K,N,nbrOfMonteCarloRealizations))/sqrt(2);
	 H = repmat(sqrt(channelVariances)',[1 N]) .* Hall(:,:,1);
	 N = size(H,2); %Number of transmit antennas (in total)


%% 遗传算法参数初始化
maxgen=20;                         %进化代数，即迭代次数
sizepop=120;                        %种群规模，即本种群中的染色体数量
T = 4;								%pool of 优秀染色体
pcross=[0.5];                       %交叉概率选择，0和1之间
pmutation=[0.09];                    %变异概率选择，0和1之间

%节点总数
numsum=N;

lenchrom=ones(1,numsum);   %初始化染色体中的基因序列长度length_of_chrom 
%需要修改bound    
bound=[-1*ones(numsum,1) 1*ones(numsum,1)];    %数据范围

%------------------------------------------------------种群初始化--------------------------------------------------------
individuals=struct('fitness',zeros(1,sizepop), 'chrom',[]);  %将种群信息定义为一个结构体
avgfitness=[];                      %每一代种群的平均适应度
bestfitness=[];                     %每一代种群的最佳适应度
bestchrom=[];                       %适应度最好的染色体
%%
individualsT=struct('fitness',zeros(1,T), 'chrom',[]);  %将挑选出来的T个种群信息定义为一个结构体
individualsTout=struct('fitness',zeros(1,sizepop-T), 'chrom',[]);  %将挑选出来的(Q-T)个种群信息定义为一个结构体
%初始化种群
for i=1:sizepop
    %随机产生一个种群
    individuals.chrom(i,:)=Code(lenchrom,bound);    %编码为实数！！！
    x=sort01(individuals.chrom(i,:));
    %计算适应度
	%需要修改fun函数--已修改
    individuals.fitness(i)=fun(x,H);   %染色体的适应度
end
FitRecord=[];
%找最好的染色体
%需要修改最好染色体的寻找方法
[bestfitness bestindex]=max(individuals.fitness);
bestchrom=individuals.chrom(bestindex,:);  %最好的染色体
avgfitness=sum(individuals.fitness)/sizepop; %染色体的平均适应度
trace=[avgfitness bestfitness]; 

%%%%
%%%%
%%%%
%%%%
for i=1:maxgen
    disp(num2str(i));%展示遗传代数
    % 选择 --修改完成版本1
    individuals=Select(individuals,sizepop); 
    avgfitness=sum(individuals.fitness)/sizepop;
    %交叉 交叉操作仅仅针对pool中的染色体
    individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound);
	sizepop = size(individuals.chrom,1);
	 for i = 1:sizepop
		individuals.chrom(i,:) = sort01(individuals.chrom(i,:));
    end
    % 变异
    individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,i,maxgen,bound);
	%sizepop = size(individuals.chrom,1);
    for i = 1:sizepop
		individuals.chrom(i,:) = sort01(individuals.chrom(i,:));
    end
    % 计算适应度 
    for j=1:sizepop
        x=individuals.chrom(j,:); %解码
        individuals.fitness(j)=fun(x,H);   
    end
    
  %找到最小和最大适应度的染色体及它们在种群中的位置
  %也需要重新设计
  
    [newbestfitness,newbestindex]=max(individuals.fitness);
	
    [worestfitness,worestindex]=min(individuals.fitness);
    % 代替上一次进化中最好的染色体
    if newbestfitness>bestfitness
        bestfitness=newbestfitness;
        bestchrom=individuals.chrom(newbestindex,:);
    end
    individuals.chrom(worestindex,:)=bestchrom;
    individuals.fitness(worestindex)=bestfitness;
    
    avgfitness=sum(individuals.fitness)/sizepop;
    
    trace=[trace;avgfitness bestfitness]; %记录每一代进化中最好的适应度和平均适应度
    %FitRecord=[FitRecord;individuals.fitness];
end
D = sort01(bestchrom);
HG =Htrans(H,D);
save('channel_GA.mat','HG');