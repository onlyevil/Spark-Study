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


修改者：韩志超和邵俊明。   	时间：2016.12.31。 
我们在“attractor.c”中实现了我们提出的社区检测算法（Attractor）。

使用步骤如下：
	（1）在第8行中更改全局变量值（我们论文中的内聚参数λ）;
	（2）在第17行中更改全局变量值（网络节点数量）;
	（3）编译代码；比如：
		在Mac：  gcc attractor.c -o attractor;
 		在linux：gcc -std=c99 attractor.c -lm -o attractor
	（4）用三个输入和一个输出执行：
		首先输入的是网络文件（其格式如下），
		第二个输入是结果文件（其格式如下），
		最后一个输入是边的数量；
	     比如：
		在bash： ./attractor Network_karate.txt Result_karate.txt 78 IMPORTANT NOTES:
	     (1)输入的网络文件的格式（见 ./Network_example.txt）：
	    	   1）每一行，代表一条无向的边，包括两个节点：源节点和目标节点
		   2）每条边只保存一次；
		   3）节点的数量从1开始，并且必须是连续的，意味着，如果网络由34个节点，最小的节点数是1，最大的节点数是34
	    （2）输出结果的格式（见 ./Result_example.txt）：每一行包括两个数字：节点数和其相应的集群
	    （3）该源代码假定最大度节点的度小于1000。如果不是这种情况，第12行的代码应该根据实际情况进行修改。

错误修复-----------sortFun函数实际上没有更新“LastChangeIndex”。
已经修复。original_table_3.jpg是原始输出。
	  correct_table_3.jpg是正确输出。