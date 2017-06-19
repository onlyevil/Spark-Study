Modified by Zhichao, Han and Junming, Shao. 2016.12.31.

We implement our proposed community detection algorithm (Attractor) in 'attractor.c'. 



Here are the steps to use it:

	(1) change the global variable value(it is the cohesion parameter $\lambda$ in our paper) in line 8;

	(2) change the global variable value(it is the network's node number) in line 17;

	(3) Compile the codes;For example: 
		On Mac: gcc attractor.c -o attractor;
			
 		On linux: gcc -std=c99 attractor.c -lm -o attractor

	(4) execute it with three input and one output: 
		first input is the network file(its format is explained below), 
		second input is the result file(its format is explained below), 
		last input is the edge number;

	    For example: 
 		In bash: ./attractor Network_karate.txt Result_karate.txt 78 


IMPORTANT NOTES:

	    (1) Format of the input network file(see ./Network_example.txt):
	
		  1) each line, represents an undirected edge, contains two numbers: source node and the target node
	
		  2) each edge is just saved only once;
	
		  3) the numbers of node begin at 1 and they must be continuous, 
		     that is, if the network has 34 nodes, the 
smallest node number is 1 and the biggest node number is 34


	    (2) Format of the output result(see ./Result_example.txt): each line contains two number: the node number and its corresponding cluster


            (3) This source code supposes that the degree of the biggest degree node is smaller than 1000. 
    		If it is not this situation, the codes in line 12 should be modified according the real situation.




Bug_Fixed
-------------------
The sortFun function did not update the "LastChangeIndex" actually. 
Already fixed.
original_table_3.jpg is the original output. 

correct_table_3.jpg is the corrected output.
-------------------


�޸��ߣ���־�����ۿ�����   	ʱ�䣺2016.12.31�� 
�����ڡ�attractor.c����ʵ���������������������㷨��Attractor����

ʹ�ò������£�
	��1���ڵ�8���и���ȫ�ֱ���ֵ�����������е��ھ۲����ˣ�;
	��2���ڵ�17���и���ȫ�ֱ���ֵ������ڵ�������;
	��3��������룻���磺
		��Mac��  gcc attractor.c -o attractor;
 		��linux��gcc -std=c99 attractor.c -lm -o attractor
	��4�������������һ�����ִ�У�
		����������������ļ������ʽ���£���
		�ڶ��������ǽ���ļ������ʽ���£���
		���һ�������Ǳߵ�������
	     ���磺
		��bash�� ./attractor Network_karate.txt Result_karate.txt 78 IMPORTANT NOTES:
	     (1)����������ļ��ĸ�ʽ���� ./Network_example.txt����
	    	   1��ÿһ�У�����һ������ıߣ����������ڵ㣺Դ�ڵ��Ŀ��ڵ�
		   2��ÿ����ֻ����һ�Σ�
		   3���ڵ��������1��ʼ�����ұ����������ģ���ζ�ţ����������34���ڵ㣬��С�Ľڵ�����1�����Ľڵ�����34
	    ��2���������ĸ�ʽ���� ./Result_example.txt����ÿһ�а����������֣��ڵ���������Ӧ�ļ�Ⱥ
	    ��3����Դ����ٶ����Ƚڵ�Ķ�С��1000��������������������12�еĴ���Ӧ�ø���ʵ����������޸ġ�

�����޸�-----------sortFun����ʵ����û�и��¡�LastChangeIndex����
�Ѿ��޸���original_table_3.jpg��ԭʼ�����
	  correct_table_3.jpg����ȷ�����