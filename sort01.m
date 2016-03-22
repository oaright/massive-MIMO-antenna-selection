function xs = sort01(x)
%input 随机产生的变量x
%output 排序好的变量xs
%代码目标 将x中数值最高的一半置1，余下的一半置0

x_length = length(x);
xs = ones(1,x_length);
[B,I] = sort(x,'descend');%I部分存的为序号
for i = 1:x_length/2
	xs(I(i)) = 1;
end
for i = x_length/2+1:x_length
	xs(I(i)) = 0;
end