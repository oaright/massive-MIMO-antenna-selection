function ret=select(individuals,sizepop，T)
% 本函数对每一代种群中的染色体进行选择，以进行后面的交叉和变异
% individuals input  : 种群信息
% sizepop     input  : 种群规模
% T           input  : pool的大小     
% ret         output : 经过选择后的种群


%%初始化
individualsT=struct('fitness',zeros(1,T), 'chrom',[]);  %将挑选出来的T个种群信息定义为一个结构体
%根据个体适应度值进行排序
[B,I] = sort(individuals.fitness,'descend');
for i = 1:T
	individualsT.chrom(i) = individuals.chrom(I(i));
	individualsT.fitness(i) = individuals.fitness(I(i));
end
ret=individualsT;
