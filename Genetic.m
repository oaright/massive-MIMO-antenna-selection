function HG = Genetic(H)
%
%input H = N x Nr 维信道矩阵
%output  HG 经过选择后的信道矩阵 

%继承参数
N = size(H,2); %Number of transmit antennas (in total)


%% 遗传算法参数初始化
maxgen=20;                         %进化代数，即迭代次数
sizepop=10;                        %种群规模，即本种群中的染色体数量
pcross=[0.2];                       %交叉概率选择，0和1之间
pmutation=[0.2];                    %变异概率选择，0和1之间

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
%初始化种群
for i=1:sizepop
    %随机产生一个种群
    individuals.chrom(i,:)=sort01(Code(lenchrom,bound));    %编码（binary和grey的编码结果为一个实数，float的编码结果为一个实数向量）
    x=individuals.chrom(i,:);
    %计算适应度
	%需要修改fun函数
    individuals.fitness(i)=fun(x,H);   %染色体的适应度
end
FitRecord=[];
%找最好的染色体
%需要修改最好染色体的寻找方法
[bestfitness bestindex]=max(individuals.fitness);
bestchrom=individuals.chrom(bestindex,:);  %最好的染色体
avgfitness=sum(individuals.fitness)/sizepop; %染色体的平均适应度
% 记录每一代进化中最好的适应度和平均适应度
trace=[avgfitness bestfitness]; 
 
%% 迭代求解最佳初始阀值和权值
% 进化开始
for i=1:maxgen
    disp(num2str(i));%展示遗传代数
    % 选择
    individuals=Select(individuals,sizepop); 
    avgfitness=sum(individuals.fitness)/sizepop;
    %交叉
    individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound);
	 for i = 1:sizepop
		individuals.chrom(i,:) = sort01(individuals.chrom(i,:));
    end
    % 变异
    individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,i,maxgen,bound);
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
    FitRecord=[FitRecord;individuals.fitness];
end
D = bestchrom;
HG =Htrans(H,D);
