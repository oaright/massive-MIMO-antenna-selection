function HG = Htrans(H,D)
%%函数功能说明：完成信道转换
%%根据选择向量D将已知信道H转换成要求的信道HG
%input: H
%input: D
%output: HG
len = length(D);
count = 1;
for i = 1: len
	if D(i) == 1
	HG(:,count) = H(:,i);
	count = count +1;
	end
end