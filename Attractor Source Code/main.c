#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// the Cohesion parameter $\lambda$ in our paper 【内聚参数】
#define BETA 0.6

// Network's node number 【网络节点数量】
#define NetNodeNum 34

// suppose max neighbors number(the biggest degree) is not bigger than 1000 【假设最大邻居数（最大度）不超过1000】
#define MaxNeighborNum 1000

// edge information
typedef struct Edge{
	int n_1_ID;	//【节点1的ID】
	int n_2_ID;	//【节点2的ID】
	int e_ID;	//【边e的ID】
	double dist;	//【距离】

	int * CN_Loc;// CN_Loc points to a 2D matrix with 4 columns 【指针，指向一个四列的二维矩阵】
	int CN_Loc_Length;// row number of CN_Loc, for search CN_Loc in Interaction 【CN_Loc的行数，用于搜索互动中的CN_Loc】

	int * N_1_Info;	//【节点1的独有邻居信息】
	int N_1_Info_Length;// row number of N_1_Info, for search N_1_Info in Interaction 【N_1_Info的行数，用于搜索互动中的N_1_Info】

	int * N_1_Loc; //【节点1的独有邻居地址】
	int N_1_Loc_Length;// row number of N_1_Loc, for search N_1_Loc in Interaction  【N_1_Loc的行数，用于搜索互动中的N_1_Loc】

	int * N_2_Info;
	int N_2_Info_Length;// row number of N_2_Info, for search N_2_Info in Interaction

	int * N_2_Loc;
	int N_2_Loc_Length;// row number of N_2_Loc, for search N_2_Loc in Interaction

}EA;

// for EdgesOfNodeTable  【节点表的边】
typedef struct EdgeLoc{
	int Row;
	int Column;
}EdgeLoc;

//【节点表头】
typedef struct PointerToNeighbors{
	int NodeNum;// 【节点的邻居数】
	int * Pointer;//【指向邻接表】
}NodeTableHead;

// (EdgeWithAttachmentTable header)【带附表头的边？】
typedef struct PointerToEdges{
	int EdgeNum;
	EA * Pointer;
}EdgeTableHead;

// (EdgesOfNodeTable header)【带节点表头的边？】
typedef struct PointerToEdgeLoc{
	int EdgeNum;
	EdgeLoc * Pointer;
}EdgesOfNodeTableHead;

/* global variable 【全局变量】*/
NodeTableHead NodeTable[NetNodeNum]; // save the network, suppose node id begin at 1【保存网络，假定节点id从1开始】
EdgeTableHead EdgeWithAttachmentTable[NetNodeNum]; // save edges with its attachments【保存有连接的边】
EdgesOfNodeTableHead EdgesOfNodeTable[NetNodeNum]; // save each node's edges(every edges that connected to it)【保存每个节点的边（与它连接的每一条边）】
int CN[MaxNeighborNum];// most neighbors number is 66【大多数邻居数是66】共同邻居
int DiffA[MaxNeighborNum];//A的独有邻居
int DiffB[MaxNeighborNum];//B的独有邻居
int NetEdgeNum;
// other global variables: EdgesOfNodeTable, clusters【其它的全局变量：节点表的边，簇】

/* Establish NodeTable【建立节点表】 */
void EstablishNodeTable(char * inputfilename, int FileLine){
    //文件名为Network_karate.txt，文件行为78，即边的数量
	FILE * file = fopen(inputfilename,"r");//将文件内容读取到file文件
	if(file==NULL)
		{printf("cannot open %s\n", inputfilename);exit(1);}
	else
		printf("opened %s\n successfully", inputfilename);

	int NeighborNumber[NetNodeNum];//邻居节点数量
	for(int i=0;i<NetNodeNum;i++)
		NeighborNumber[i] = 0;//对于每个节点i，初始邻居节点数量为0
	int node_1;//用于保存文件中每条边的节点
	int node_2;

	for(int i=0;i<FileLine;i++){//对于文件的每一行，即每条边
		fscanf(file,"%d %d",&node_1,&node_2);//读取边的两个节点，保存到node_1,node_2
		NeighborNumber[node_1-1]++;//有边则该节点的邻居节点+1
		NeighborNumber[node_2-1]++;
	}//用于求每个节点的邻居节点数

	//initiate NodeTable【初始化节点表】
	for(int i=0;i<NetNodeNum;i++){//对于每个网络节点
		NodeTable[i].NodeNum = 0; //NodeTable属于结构体NodeTableHead,初始化每个节点的邻居数为0
		NodeTable[i].Pointer = (int *)malloc(NeighborNumber[i]*sizeof(int));
		//为每个节点的邻居数组申请空间
		if(NodeTable[i].Pointer==NULL)
			{printf("Memory is not enough when initiating NodeTable\n");exit(1);}
	}
	fclose(file);

	file = fopen(inputfilename,"r");//再次打开输入文件，从文件首行重新读取
	if(file==NULL)
		{printf("cannot open %s\n", inputfilename);exit(1);}
	else
		printf("have opened %s\n", inputfilename);

	for(int i=0;i<FileLine;i++){
		fscanf(file,"%d %d",&node_1,&node_2);
		//node_1和node_2为一条边，node_2为node_1的邻居节点
		NodeTable[node_1-1].Pointer[NodeTable[node_1-1].NodeNum++] = node_2;
		//在节点和其邻居节点间建立连接
		NodeTable[node_2-1].Pointer[NodeTable[node_2-1].NodeNum++] = node_1;
		//因为文件中节点id是从1开始，而代码中是从0开始，所以需要-1
	}//用于建立每个节点的节点表
	// check whether NodeTable[i].NodeNum equals to NeighborNumber[i] or not【检查NodeTable[i].NodeNum是否等于NeighborNumber[i]】
	for(int i=0;i<NetNodeNum;i++)
		if(NodeTable[i].NodeNum!=NeighborNumber[i])//用于判断节点表是否建立错误
			{printf("NodeTable[%d] error\n", i+1);exit(1);}
	fclose(file);


	// sort each node's neighbors in ascending order【按照升序为每个节点的邻居排序】
	void SortFun(int * , int );//declaration【声明】
	for(int i=0;i<NetNodeNum;i++)
		SortFun(NodeTable[i].Pointer, NodeTable[i].NodeNum);

}

/* Establish EdgeTable 【建立边表】*/
void EstablishEdgeTable(){
	int EdgeNumber[NetNodeNum];// save edge number of each node【保存每个节点的边数】
	//compute node[loop_i]'s neighbors number【计算node[loop_i]的邻居数】
	for(int loop_i=0;loop_i<NetNodeNum;loop_i++){//对于每个网络中的每个节点
		EdgeNumber[loop_i] = 0;//初始化每个节点的边数为0
		for(int loop_j=0; loop_j<NodeTable[loop_i].NodeNum; loop_j++)//对于每个节点的邻居节点
			if(NodeTable[loop_i].Pointer[loop_j]>loop_i+1)// save each edge only once【每个边只保存一次】
			//觉得应该是>=吧   【之前建立节点表时是边重复建了一次】
			//明白了，确实是>,因为pointer指向的节点（从1开始),而loop_i是从0开始的
				EdgeNumber[loop_i]++;
	}

	// compute n_1_ID, n_2_ID, e_ID
	int edgeid = 0;
	for(int loop_i=0;loop_i<NetNodeNum;loop_i++){   //对于网络中每个节点
		EdgeWithAttachmentTable[loop_i].EdgeNum = EdgeNumber[loop_i];//每个节点的边的数量
		//EdgeWithAttachmentTable 为结构体 EdgeTableHead
		EdgeWithAttachmentTable[loop_i].Pointer = (EA *)malloc(EdgeNumber[loop_i]*sizeof(EA));
		//为每个节点申请【边的数量*EA】的空间，EA为结构体（Pointer占用该结构体的空间）
		if(EdgeWithAttachmentTable[loop_i].Pointer==NULL)
			{printf("memory is not enough with node%d's edges\n", loop_i+1);exit(1);}

		int edge_loc_tmp=0;//边表中边的下标（从0开始）
		for(int loop_j=0; loop_j<NodeTable[loop_i].NodeNum; loop_j++){
        //对于网络中每一个节点loop_i的邻居节点loop_j
			if(NodeTable[loop_i].Pointer[loop_j]>loop_i+1){//同上，去掉重复边
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].n_1_ID = loop_i+1;
                //下标为loop_i(下标从0开始）的节点的第edge_loc_tmp条边的一个节点为loop_i+1（节点从1开始）
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].n_2_ID = NodeTable[loop_i].Pointer[loop_j];
				//下标为loop_i(下标从0开始）的节点的第edge_loc_tmp条边的另一个节点从loop_i的节点表中获得
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].e_ID = edgeid;
				//第loop_i个节点的第edge_loc_tmp条边的id为 edgeid（从0开始）
				edgeid++;//表示的是网络中总的边数
				NetEdgeNum = edgeid;
				edge_loc_tmp++;
			}
		}
		// check whether edge_loc_tmp equals to EdgeWithAttachmentTable[loop_i].EdgeNum or not【检查edge_loc_tmp是否等于dgeWithAttachmentTable[loop_i].EdgeNum】
		if(edge_loc_tmp!=EdgeWithAttachmentTable[loop_i].EdgeNum)
			{printf("EdgeWithAttachmentTable[%d] error(line 175)\n", loop_i+1); exit(1);}
	}

	/* compute dist, CN_Loc, N_1_Info, N_1_Loc, N_2_Info, N_2_Loc */
	for(int i=0;i<NetNodeNum;i++){
		for(int j=0;j<EdgeWithAttachmentTable[i].EdgeNum;j++){
			int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
			int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
			//node_1和node_2分别代表第i个节点第j条边的两个节点
			void FindNeighbors(int, int);// function declaration【函数声明】
			FindNeighbors(node_1,node_2);//发现两个节点的共同邻居，以及各自的独有邻居


			// compute distance【计算距离】
			//第i个节点的第j条边的杰卡德距离，CN[0]为共同邻域的长度、diffA[0]、diffB[0]也是总长度
			EdgeWithAttachmentTable[i].Pointer[j].dist = 1.0 - (double)(CN[0]+2)/(double)(CN[0]+DiffA[0]+DiffB[0]+2);


			/* establish CN_Loc 【四列的二维矩阵】建来做什么？？？用于保存每个公共邻居点的四列信息*/
			//第一列：node_1和nodecn中最小节点对应id，行
			//第二列：node_1和NodeCN所在边对应的id（相对于最小节点），列
			//第三列：node_2和nodecn中最小节点对应id
			//第四列：node_2和NodeCN所在边对应的id（相对于最小节点）
			EdgeWithAttachmentTable[i].Pointer[j].CN_Loc = (int *)malloc(4*CN[0]*sizeof(int));
			if(EdgeWithAttachmentTable[i].Pointer[j].CN_Loc==NULL)
				{printf("EdgeWith AttachmentTable[%d].Pointer is NULL\n", i);exit(1);}
			EdgeWithAttachmentTable[i].Pointer[j].CN_Loc_Length = CN[0];

			for(int s=1;s<=CN[0];s++){//对于第i个节点的第j条边的共同邻居点
				int NodeCN = CN[s];//NodeCN表示第s个共同邻居点
				int Loc_Node_1_R = -1;//node_1和nodecn中最小节点对应id
				int Loc_Node_1_C = -1;//node_1和NodeCN所在边对应的id（相对于最小节点）
				int Loc_Node_2_R = -1;//node_2和nodecn中最小节点对应id
				int Loc_Node_2_C = -1;//node_2和NodeCN所在边对应的id（相对于最小节点）
				//找出node_1、NodeCN中的最大节点和最小节点
				int NodeMin = (node_1<NodeCN)?node_1:NodeCN;
				int NodeMax = (node_1>NodeCN)?node_1:NodeCN;
				Loc_Node_1_R = NodeMin - 1;//节点id转为下标值要注意-1。表示的是最小节点的下标值

				//【找出node_1--NodeCN 边的id，用Loc_Node_1_C表示】
				for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
                    //对于最小节点的第loop条边
                    //如果第loop条边的n_2_ID恰等于最大节点
					if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
						{	Loc_Node_1_C = loop;	break;}

				if(Loc_Node_1_C==-1)
					{printf("line 207 error\n"); exit(1);}

				NodeMin = (node_2<NodeCN)? node_2:NodeCN;
				NodeMax = (node_2>NodeCN)? node_2:NodeCN;
				Loc_Node_2_R = NodeMin - 1;
                //同上
				for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
					if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
						{	Loc_Node_2_C = loop;	break;}
				if(Loc_Node_2_C==-1)
					{printf("line 261 error\n"); exit(1);}
                //s为公共邻居点的下标，从1开始【针对每一个公共邻居点，占据4行1列】
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4*(s-1)+0] = Loc_Node_1_R;
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4*(s-1)+1] = Loc_Node_1_C;
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4*(s-1)+2] = Loc_Node_2_R;
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4*(s-1)+3] = Loc_Node_2_C;
			}

//			printf("Node %d, Node %d complete CN_Loc\n", node_1, node_2);



			/* establish N_1_Loc, N_1_Info【建立Node_1的独有邻居】*/
			//【N_1_Info】
			//Node_N_1，节点1的独有邻居
			//CN[0]，节点1的独有邻居和节点2的共同邻居长度
			//节点1和节点1的独有邻居中最小的节点的下标，也即在边表中的行
			//初始为-1，表示节点1和节点1的独有邻居所在边的id（根据最小节点），也即在边表中的列
			EdgeWithAttachmentTable[i].Pointer[j].N_1_Info = (int *)malloc(2*2*DiffA[0]*sizeof(int));
			//【第i个节点的第j条边的节点1的独有邻居的信息】
			// each node have four elements: node id, the number common neighbors of N_1 and node_2, row of <N_1, node_1>, column of <N_1, node_1>
			//【每个节点有四个元素：节点id，N_1和node_2的共同邻域数，<N_1,node_1>的行，<N_1,node_1>的列
			EdgeWithAttachmentTable[i].Pointer[j].N_1_Info_Length = DiffA[0];
            //【第i个节点的第j条边的节点1的独有邻居的长度为DiffA[0]】

			if(EdgeWithAttachmentTable[i].Pointer[j].N_1_Info==NULL)
				{printf("memory not enough With N_1_Info\n"); exit(1);}

			int Edgenum_tmp = 0;//
			for(int s=1;s<=DiffA[0];s++){//对于节点1的每一个独有邻居
				int Node_N_1 = DiffA[s];//Node_N_1表示节点1的第s个独有邻居
				void FindCN(int, int);// function declaration函数说明
				FindCN(Node_N_1, node_2);// CN changes【CN改变，独有邻居和节点2的公共邻居】
				Edgenum_tmp  = Edgenum_tmp + CN[0];//节点1的所有独有邻居与节点2的公共邻居长度
			}

			//节点1独有邻居与节点2的共同邻居的信息表N_1_Loc
			//最小节点下标，也即<Node_N_1,NodeCN>所在行(边表)
			//<Node_N_1,NodeCN>所在列(边表)
			//<node_2,NodeCN>所在行(边表)
			//<node_2,NodeCN>所在列(边表)
			EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc = (int *)malloc(4*Edgenum_tmp*sizeof(int));
			if(EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc==NULL)
				{printf("memory not enough With N_1_Loc(line 243)\n"); exit(1);}
			EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc_Length = Edgenum_tmp;

			Edgenum_tmp = 0;
			for(int s=1;s<=DiffA[0];s++){
				int Node_N_1 = DiffA[s];
				void FindCN(int,int);
				FindCN(Node_N_1, node_2);//CN changes

				// N_1_Info[4*(s-1)+0] is N_1 id, N_1_Info[4*(s-1)+1] is the number of common neighbors of N_1 and Node_2
				//【N_1_Info[4*(s-1)+0]是N_1 id，N_1_Info[4*(s-1)+1]是N_1和Node_2的共同邻域数】
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+0] = Node_N_1;// this node's id
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+1] = CN[0];// number or common neighbors of this node and node_2

				int NodeMin;
				int NodeMax;

				// N_1_Info[4*(s-1)+2] is the location of row of <N_1, node_1>, N_1_Info[4*(s-1)+3] is the column of <N_1, node_1>
				//【N_1_Info[4*(s-1)+2]】是<N_1, node_1>的行的位置，N_1_Info[4*(s-1)+3]是<N_1, node_1>的列
				NodeMin = (node_1<Node_N_1)?node_1:Node_N_1;
				NodeMax = (node_1>Node_N_1)?node_1:Node_N_1;
				//节点1和节点1的独有邻居大小比较

				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+2] = NodeMin-1;//row
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+3] = -1;

				//找出节点1和节点1的独有邻居所在边的id（根据最小节点）
				for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
					if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
						{	EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+3] = loop;	break;}
				if(EdgeWithAttachmentTable[NodeMin-1].EdgeNum==-1)
					{printf("line 270 error\n");exit(1);}

				for(int ss=1;ss<=CN[0];ss++){//CN[0]表示节点1独有邻居和节点2的共同邻居长度
					int NodeCN = CN[ss];//NodeCN表示共同邻居节点
					int Loc_Node_N_1_R = -1;//最小节点下标，也即<Node_N_1,NodeCN>所在行(边表)
					int Loc_Node_N_1_C = -1;//<Node_N_1,NodeCN>所在列(边表)
					int Loc_Node_2_R = -1;//<node_2,NodeCN>所在行(边表)
					int Loc_Node_2_C = -1;//<node_2,NodeCN>所在列(边表)

					//节点1的独有邻居和共同邻居节点比较大小
					NodeMin = (Node_N_1<NodeCN)?Node_N_1:NodeCN;
					NodeMax = (Node_N_1>NodeCN)?Node_N_1:NodeCN;

					Loc_Node_N_1_R = NodeMin - 1;
					for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
						if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
							{	Loc_Node_N_1_C = loop;	break;}

					NodeMin = (node_2<NodeCN)?node_2:NodeCN;
					NodeMax = (node_2>NodeCN)?node_2:NodeCN;
					Loc_Node_2_R = NodeMin - 1;
					for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
						if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
							{	Loc_Node_2_C = loop;	break;}

					EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc[4*Edgenum_tmp+0] = Loc_Node_N_1_R;
					EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc[4*Edgenum_tmp+1] = Loc_Node_N_1_C;
					EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc[4*Edgenum_tmp+2] = Loc_Node_2_R;
					EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc[4*Edgenum_tmp+3] = Loc_Node_2_C;
					Edgenum_tmp++;
				}
			}


			/* establish N_2_Loc, N_2_Info 【建立N_2_Loc, N_2_Info】*/
			EdgeWithAttachmentTable[i].Pointer[j].N_2_Info = (int *)malloc(2*2*DiffB[0]*sizeof(int));
			EdgeWithAttachmentTable[i].Pointer[j].N_2_Info_Length = DiffB[0];

			Edgenum_tmp = 0;
			for(int s=1;s<=DiffB[0];s++){
				int Node_N_2 = DiffB[s];
				void FindCN(int,int);// declaration
				FindCN(Node_N_2, node_1);//CN changes
				Edgenum_tmp = Edgenum_tmp + CN[0];
			}
			EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc = (int *)malloc(4*Edgenum_tmp*sizeof(int));
			EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc_Length = Edgenum_tmp;

			Edgenum_tmp = 0;
			for(int s=1;s<=DiffB[0];s++){
				int Node_N_2 = DiffB[s];
				void FindCN(int, int);
				FindCN(Node_N_2, node_1);// CN changes

				// N_2_Info
				EdgeWithAttachmentTable[i].Pointer[j].N_2_Info[4*(s-1)+0] = Node_N_2;
				EdgeWithAttachmentTable[i].Pointer[j].N_2_Info[4*(s-1)+1] = CN[0];

				int NodeMin;
				int NodeMax;

				NodeMin = (node_2<Node_N_2)?node_2:Node_N_2;
				NodeMax = (node_2>Node_N_2)?node_2:Node_N_2;
				EdgeWithAttachmentTable[i].Pointer[j].N_2_Info[4*(s-1)+2] = NodeMin-1;
				for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
					if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
						{	EdgeWithAttachmentTable[i].Pointer[j].N_2_Info[4*(s-1)+3] = loop;	break;}

				for(int ss=1;ss<=CN[0];ss++){
					int NodeCN = CN[ss];
					int Loc_Node_N_2_R;
					int Loc_Node_N_2_C;
					int Loc_Node_1_R;
					int Loc_Node_1_C;
					NodeMin = (Node_N_2<NodeCN)?Node_N_2:NodeCN;
					NodeMax = (Node_N_2>NodeCN)?Node_N_2:NodeCN;
					Loc_Node_N_2_R = NodeMin - 1;
					for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
						if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
							{	Loc_Node_N_2_C = loop;	break;}

					NodeMin = (node_1<NodeCN)?node_1:NodeCN;
					NodeMax = (node_1>NodeCN)?node_1:NodeCN;
					Loc_Node_1_R = NodeMin-1;
					for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
						if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
							{	Loc_Node_1_C = loop;	break;}

					EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc[4*Edgenum_tmp+0] = Loc_Node_N_2_R;
					EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc[4*Edgenum_tmp+1] = Loc_Node_N_2_C;
					EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc[4*Edgenum_tmp+2] = Loc_Node_1_R;
					EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc[4*Edgenum_tmp+3] = Loc_Node_1_C;
					Edgenum_tmp++;
				}
				// check Edgenum_tmp equals to EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc's length
				//【检查Edgenum_tmp是否等于EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc的长度】
			}
		}
	}

	// FILE * OUT = fopen("result/edgenum.txt","w");
	// for(int i=0;i<NetNodeNum;i++){
	// 	fprintf(OUT,"Node %d: Edge Neighbors number:%d\n", i+1, EdgeWithAttachmentTable[i].EdgeNum);
	// }
	// fclose(OUT);
}

void EstablishEdgesOfNodeTable(){
	for(int i=0;i<NetNodeNum;i++){
		EdgesOfNodeTable[i].EdgeNum = 0;//结构体 PointerToEdgeLoc
		EdgesOfNodeTable[i].Pointer = (EdgeLoc *)malloc(NodeTable[i].NodeNum*sizeof(EdgeLoc));
		if(EdgesOfNodeTable[i].Pointer==NULL)
			{printf("memory not enough with EdgesOfNodeTable(line 378)\n");exit(1);}
	}

	for(int i=0;i<NetNodeNum;i++)
		for(int j=0;j<EdgeWithAttachmentTable[i].EdgeNum;j++){//属于结构体 PointerToEdges EdgeTableHead，每个节点的边的数量
			int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
			int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
			// a little check
			if(node_1!=i+1)
				{printf("line 386 error\n"); exit(1);}

            //第i个节点的第EdgeNum个节点所在边的行
			EdgesOfNodeTable[node_1-1].Pointer[EdgesOfNodeTable[node_1-1].EdgeNum].Row = i;
			EdgesOfNodeTable[node_1-1].Pointer[EdgesOfNodeTable[node_1-1].EdgeNum].Column = j;
			EdgesOfNodeTable[node_1-1].EdgeNum++;

			EdgesOfNodeTable[node_2-1].Pointer[EdgesOfNodeTable[node_2-1].EdgeNum].Row = i;
			EdgesOfNodeTable[node_2-1].Pointer[EdgesOfNodeTable[node_2-1].EdgeNum].Column = j;
			EdgesOfNodeTable[node_2-1].EdgeNum++;
		}

	for(int i=0;i<NetNodeNum;i++)
		if(EdgesOfNodeTable[i].EdgeNum!=NodeTable[i].NodeNum)
			{printf("node_%d: Edges of this node does not math neighbors of this node(line 398 error)\n", i+1); exit(1);}
		else
			printf("Node: %d neighbors number: %d\n", i+1, EdgesOfNodeTable[i].EdgeNum);
}

/* Interaction 【交互】*/
void Interaction(double beta, int FileLine){
	// beta controls clusters number【β控制集群数量】，FileLine为边的数量，手动输入

	if(NetEdgeNum!=FileLine){
		printf("NetEdgeNum error...(line 413)\n");
		exit(1);
	}
	double D[NetEdgeNum];//边的距离

	// esablish D
	int EdgeLocCount = 0;
	for(int s_1=0;s_1<NetNodeNum;s_1++)
		for(int s_2=0;s_2<EdgeWithAttachmentTable[s_1].EdgeNum;s_2++)
			D[EdgeLocCount++] = EdgeWithAttachmentTable[s_1].Pointer[s_2].dist;
	if(EdgeLocCount!=NetEdgeNum)
		{printf("line 424 error\n");exit(1);}
    // for(int s_1=0;s_1<NetNodeNum;s_1++)
    //     printf("%f\n", D[s_1]);

	// Interaction process
	int Terminate = 1;//终止条件
	int Loop = 0;
	while(Terminate){
		Loop++;//迭代次数
		printf("Loop: %d\n", Loop);//第一次迭代
		for(int s_1=0;s_1<NetNodeNum;s_1++)
			for(int s_2=0;s_2<EdgeWithAttachmentTable[s_1].EdgeNum;s_2++)
				if(EdgeWithAttachmentTable[s_1].Pointer[s_2].dist>0 && EdgeWithAttachmentTable[s_1].Pointer[s_2].dist<1){
                //如果第s_1个节点的第s_2条边的距离（0,1）
					EA ThisEA = EdgeWithAttachmentTable[s_1].Pointer[s_2];
                    //用ThisEA表示这条边
					int ThisNode_1 = ThisEA.n_1_ID;
					int ThisNode_2 = ThisEA.n_2_ID;
                    int ThisEdgeId = ThisEA.e_ID;

					double CI = 0.0;// common neighbors' influence【共同邻域的影响】
					double N_1_I = 0.0;// node_1's neighbors' influence【node_1的独有邻居影响】
					double N_2_I = 0.0;// node_2's neighbors' influence【node_2的独有邻居影响】

					/* common neighbors' influence【共同邻域影响】 */
					for(int s_3=0;s_3<ThisEA.CN_Loc_Length;s_3++){//共同邻域的长度
						double Distance_CI_1;//节点1和共同邻居节点间的距离
						double Distance_CI_2;//节点2和共同邻居节点间的距离
						int * CNLocTmp = ThisEA.CN_Loc;//共同邻域表

						//节点1和共同邻居节点间的距离
						Distance_CI_1 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+0]].Pointer[CNLocTmp[s_3*4+1]].dist;
						Distance_CI_2 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+2]].Pointer[CNLocTmp[s_3*4+3]].dist;

						//<节点1,共同邻居节点>所在边的两个节点
						int CNTmp_1 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+0]].Pointer[CNLocTmp[s_3*4+1]].n_1_ID;
						int CNTmp_2 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+0]].Pointer[CNLocTmp[s_3*4+1]].n_2_ID;

						int CNTmp_3 = (CNTmp_1==ThisNode_1)?CNTmp_2:CNTmp_1;//表示邻居节点
						CI = CI + ( (double)(1-Distance_CI_2) / (double)NodeTable[ThisNode_1-1].NodeNum ) *  sin(1-Distance_CI_1);
						// printf("CI influence1: %f\n", ( (double)(1-Distance_CI_2) / (double)NodeTable[ThisNode_1-1].NodeNum ) *  sin(1-Distance_CI_1));


						//<节点2,共同邻居节点>所在边的两个节点
						CNTmp_1 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+2]].Pointer[CNLocTmp[s_3*4+3]].n_1_ID;
						CNTmp_2 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+2]].Pointer[CNLocTmp[s_3*4+3]].n_2_ID;

						CNTmp_3 = (CNTmp_1==ThisNode_2)?CNTmp_2:CNTmp_1;//没用到
						CI = CI + ( (double)(1-Distance_CI_1) / (double)NodeTable[ThisNode_2-1].NodeNum ) *  sin(1-Distance_CI_2);
						// printf("CI influence2: %f\n", ( (double)(1-Distance_CI_1) / (double)NodeTable[ThisNode_2-1].NodeNum ) *  sin(1-Distance_CI_2));
					}

					/* node_1's neighbors' influence 【节点1的独有邻居影响】*/
					int s_3 = 0;//节点1的独有邻居
					int s_4 = 0;
					int s_4_count;
					while(s_3<ThisEA.N_1_Info_Length){//节点1的独有邻居的长度（从0开始）
						int Ngh_N_1 = ThisEA.N_1_Info[s_3*4+0];// node_1's neighbor's id【节点1独有邻居节点】
						int CN_Num = ThisEA.N_1_Info[s_3*4+1];// common neighbors' number of node_1's neighbor and node_2
						//【<节点1独有邻居，节点2>共同邻居长度】

						//ThisEA.N_1_Info[s_3*4+2]指<节点1，节点1独有邻居>在边表中的行
						//ThisEA.N_1_Info[s_3*4+3]指<节点1，节点1独有邻居>在边表中的列
						double Distance_N_N_1 = EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3*4+2]].Pointer[ThisEA.N_1_Info[s_3*4+3]].dist;
                        //<节点1，节点1独有邻居>的距离

						s_4_count = s_4;//表示的是节点1的所有独有邻居与节点2的共同邻居表下标
						double lambda_numerator = 0;//分子
						//一个节点1的独有邻居的lambda=节点1独有邻居与所有共同邻居的lambda之和+节点2与所有共同邻居的lambda之和
						for(s_4=s_4_count;s_4<s_4_count+CN_Num;s_4++){//【<节点1独有邻居，节点2>共同邻居长度】
							int loc_1;
							int loc_2;

							loc_1 = ThisEA.N_1_Loc[s_4*4+0];//<Node_N_1,NodeCN>所在行（边表）
							loc_2 = ThisEA.N_1_Loc[s_4*4+1];//<Node_N_1,NodeCN>所在列（边表）
							lambda_numerator = lambda_numerator + (1-EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
							loc_1 = ThisEA.N_1_Loc[s_4*4+2];//<Node_2,NodeCN>所在行（边表）
							loc_2 = ThisEA.N_1_Loc[s_4*4+3];//<Node_2,NodeCN>所在列（边表）
							lambda_numerator = lambda_numerator + (1-EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
						}

						double lambda_denominator = 0;//分母=节点1的独有邻居的边的lambda+节点2的边的lambda
						for(int s_5=0;s_5<EdgesOfNodeTable[Ngh_N_1-1].EdgeNum;s_5++){//节点1独有邻居的边数
							int loc_r = EdgesOfNodeTable[Ngh_N_1-1].Pointer[s_5].Row;//节点1独有邻居的第s_5条边的行
							int loc_c = EdgesOfNodeTable[Ngh_N_1-1].Pointer[s_5].Column;//节点1独有邻居的第s_5条边的列
							lambda_denominator = lambda_denominator + (1-EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}

						for(int s_5=0;s_5<EdgesOfNodeTable[ThisNode_2-1].EdgeNum;s_5++){
							int loc_r = EdgesOfNodeTable[ThisNode_2-1].Pointer[s_5].Row;
							int loc_c = EdgesOfNodeTable[ThisNode_2-1].Pointer[s_5].Column;
							lambda_denominator = lambda_denominator + (1-EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}

						double lambda = (double)lambda_numerator/(double)lambda_denominator;
						// double max = (lambda>beta)?lambda:beta;
						double parameter = (lambda>beta)?1:-1;
						if(parameter==1)
							N_1_I = N_1_I + parameter*lambda/(double)NodeTable[ThisNode_1-1].NodeNum * sin(1 - Distance_N_N_1);
						else
							N_1_I = N_1_I + parameter*(beta-lambda)/(double)NodeTable[ThisNode_1-1].NodeNum * sin(1 - Distance_N_N_1);

						s_3++;
					}


					/* node_2's neighbors' influence */
					s_3 = 0;
					s_4 = 0;
					while(s_3<ThisEA.N_2_Info_Length){
						int Ngh_N_2 = ThisEA.N_2_Info[s_3*4+0];
						int CN_Num = ThisEA.N_2_Info[s_3*4+1];
						double Distance_N_N_2 = EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3*4+2]].Pointer[ThisEA.N_2_Info[s_3*4+3]].dist;

						s_4_count = s_4;
						double lambda_numerator = 0;
						for(s_4=s_4_count;s_4<s_4_count+CN_Num;s_4++){
							int loc_1;
							int loc_2;

							loc_1 = ThisEA.N_2_Loc[s_4*4+0];
							loc_2 = ThisEA.N_2_Loc[s_4*4+1];
							lambda_numerator = lambda_numerator + (1-EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
							loc_1 = ThisEA.N_2_Loc[s_4*4+2];
							loc_2 = ThisEA.N_2_Loc[s_4*4+3];
							lambda_numerator = lambda_numerator + (1-EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
						}

						double lambda_denominator = 0;
						for(int s_5=0;s_5<EdgesOfNodeTable[Ngh_N_2-1].EdgeNum;s_5++){
							int loc_r = EdgesOfNodeTable[Ngh_N_2-1].Pointer[s_5].Row;
							int loc_c = EdgesOfNodeTable[Ngh_N_2-1].Pointer[s_5].Column;
							lambda_denominator = lambda_denominator + (1-EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						for(int s_5=0;s_5<EdgesOfNodeTable[ThisNode_1-1].EdgeNum;s_5++){
							int loc_r = EdgesOfNodeTable[ThisNode_1-1].Pointer[s_5].Row;
							int loc_c = EdgesOfNodeTable[ThisNode_1-1].Pointer[s_5].Column;
							lambda_denominator = lambda_denominator + (1-EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}

						double lambda = (double)lambda_numerator/(double)lambda_denominator;
						// double max = (lambda>beta)?lambda:beta;
						double parameter = (lambda>beta)?1:-1;
						if(parameter==1)
							N_2_I = N_2_I + parameter *lambda / (double)NodeTable[ThisNode_2-1].NodeNum * sin(1 - Distance_N_N_2);
						else
							N_2_I = N_2_I + parameter *(beta-lambda)/ (double)NodeTable[ThisNode_2-1].NodeNum * sin(1 - Distance_N_N_2);

						s_3++;
					}

					//第loop次迭代后，三种影响模式后的所有边的距离
					D[ThisEdgeId] = D[ThisEdgeId] + (-(N_1_I + N_2_I)-( ((double)1.0/(double)NodeTable[ThisNode_1-1].NodeNum + (double)1.0/(double)NodeTable[ThisNode_2-1].NodeNum) * sin(1-ThisEA.dist) +CI));
					if(D[ThisEdgeId]<0)
						D[ThisEdgeId] = 0;
					if(D[ThisEdgeId]>1)
						D[ThisEdgeId] = 1;
				}

		double sum_1 = 0;
		double sum_2 = 0;
		int EdgeCounter = 0;
		for(int s_1=0;s_1<NetNodeNum;s_1++)
			for(int s_2=0;s_2<EdgeWithAttachmentTable[s_1].EdgeNum;s_2++){
				sum_1 = sum_1 + EdgeWithAttachmentTable[s_1].Pointer[s_2].dist;
				sum_2 = sum_2 + D[EdgeCounter];
				EdgeCounter++;
			}
		if(sum_1==sum_2 || Loop>1000)//end condition【结束条件】迭代次数超过1000次，或迭代后边的距离没有改变
			Terminate = 0;

		// renew edge's distance【更新边的距离】
		EdgeCounter = 0;
		for(int s_1=0;s_1<NetNodeNum;s_1++)
			for(int s_2=0;s_2<EdgeWithAttachmentTable[s_1].EdgeNum;s_2++){
				EdgeWithAttachmentTable[s_1].Pointer[s_2].dist = D[EdgeCounter];
				EdgeCounter++;
			}

	}

}

/* find clusters, by BFS【通过广度优先搜索算法发现簇】 */
void FindClusters(char * outputfile){
    //结果输出文件名
	int CountNeighborNumber[NetNodeNum];//节点邻居数
	// establish CountNeighborNumber【建立计数邻居号】
	for(int i=0;i<NetNodeNum;i++)
		CountNeighborNumber[i] = 0;
	for(int i=0;i<NetNodeNum;i++)
		for(int j=0;j<EdgeWithAttachmentTable[i].EdgeNum;j++)
			if(EdgeWithAttachmentTable[i].Pointer[j].dist<1){//等于1的边不要
				int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
				int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
				CountNeighborNumber[node_1-1]++;
				CountNeighborNumber[node_2-1]++;
			}
	// establish new NodeTable【建立新的节点表，距离均为0】
	NodeTableHead NodeTable_New[NetNodeNum];
	for(int i=0;i<NetNodeNum;i++){
		NodeTable_New[i].NodeNum = 0;
		NodeTable_New[i].Pointer = (int *)malloc(CountNeighborNumber[i]*sizeof(int));
		if(NodeTable_New[i].Pointer==NULL){
			printf("memory is not enough with NodeTable_New\n");
			exit(1);
		}
	}
	for(int i=0;i<NetNodeNum;i++)//建立新的节点表
		for(int j=0;j<EdgeWithAttachmentTable[i].EdgeNum;j++)
			if(EdgeWithAttachmentTable[i].Pointer[j].dist<1){
				int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
				int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
				NodeTable_New[node_1-1].Pointer[NodeTable_New[node_1-1].NodeNum] = node_2;
				NodeTable_New[node_2-1].Pointer[NodeTable_New[node_2-1].NodeNum] = node_1;
				NodeTable_New[node_1-1].NodeNum++;
				NodeTable_New[node_2-1].NodeNum++;
			}

	//check if NodeTable is ok【检查节点表是否正常】
	for(int i=0;i<NetNodeNum;i++)
		if(NodeTable_New[i].NodeNum!=CountNeighborNumber[i])
			printf("NodeTable_New Error! In Node: %d\n", i+1);

	/* find clusters */
	int clusters[NetNodeNum+1];// clusters[i] is node_i's id, clusters[0] is unused【cluster[i]是node_i的id，cluster[0]未被使用】
	for(int k=0;k<=NetNodeNum;k++)
		clusters[k] = -1;
	int Q[NetNodeNum];// the queue【队列】
	int ClusterID = 1;
	int Terminate = 1;//终止条件

	while(Terminate){
		Terminate = 0;
		int id;
		for(id=1;id<=NetNodeNum;id++)
			if(clusters[id]==-1){//开始每个节点的簇均为-1，表示该节点尚未分簇
				Terminate = 1;//有尚未分簇的节点便不能终止循环
				clusters[id] = ClusterID;
				Q[0] = id;
				int first = -1;
				int last = 0;
				int v;
				while(first!=last){
					v = Q[++first];
					for(int len=0;len<NodeTable_New[v-1].NodeNum;len++){
						int RelatedNode = NodeTable_New[v-1].Pointer[len];
						if(clusters[RelatedNode]==-1){
							Q[++last] = RelatedNode;
							clusters[RelatedNode] = ClusterID;
						}
					}
				}
				ClusterID = ClusterID + 1;
			}
	}

	// save clusters in 2D array【将簇保存在二维数组】
	FILE * fout = fopen(outputfile,"w");
	if(fout==NULL)
		{	printf("opening outputfile fails\n"); exit(0);}
	for(int i=1;i<=NetNodeNum;i++)
		fprintf(fout, "%d %d\n", i, clusters[i]);//网络节点数，所在簇
	fclose(fout);



}

/* compute cn, diffA, diffB【计 算cn，diffA，diffB】 */
void FindNeighbors(int node_i, int node_j){

	int num1 = NodeTable[node_i-1].NodeNum;//表示node_i的邻居总数
	int * A1 = NodeTable[node_i-1].Pointer;//指向node_i的邻接表
	int num2 = NodeTable[node_j-1].NodeNum;//同上
	int * A2 = NodeTable[node_j-1].Pointer;

	int p1_loc = 0;//表示node_1节点的邻接表下标
	int p2_loc = 0;
	int * p1 = &A1[0];//p1指向node_i的第一个邻接点
	int * p2 = &A2[0];

	int cn_length = 0;//邻域的大小
	int diffa_length = 0;//node_i独有邻居的的大小
	int diffb_length = 0;

	int cn_loc = 1;//邻域地址下标，从1开始，0用于保存大小长度
	int diffa_loc = 1;//node_i独有邻居的下标，从1开始
	int diffb_loc = 1;

	while(p1_loc<num1 && p2_loc<num2){
		if(A1[p1_loc]<A2[p2_loc]){
			if(A1[p1_loc]!=node_j){
				DiffA[diffa_loc] = A1[p1_loc];
				diffa_length++;
				diffa_loc++;
				p1_loc++;
			}
			else{
				p1_loc++;
			}

		}
		else if(A1[p1_loc]==A2[p2_loc]){
			CN[cn_loc] = A1[p1_loc];//CN[]下标从1开始
			cn_length++;//邻域长度+1
			cn_loc++;
			p1_loc++;
			p2_loc++;
		}
		else{
			if(A2[p2_loc]!=node_i){
				DiffB[diffb_loc] = A2[p2_loc];
				diffb_length++;
				diffb_loc++;
				p2_loc++;
			}
			else
				p2_loc++;
		}
	}
	if(p1_loc==num1){
		while(p2_loc<num2){
			if(A2[p2_loc]!=node_i){
				DiffB[diffb_loc] = A2[p2_loc];
				diffb_length++;
				diffb_loc++;
				p2_loc++;
			}
			else
				p2_loc++;
		}
	}
	else{
		while(p1_loc<num1){
			if(A1[p1_loc]!=node_j){
				DiffA[diffa_loc] = A1[p1_loc];
				diffa_length++;
				diffa_loc++;
				p1_loc++;
			}
			else
				p1_loc++;
		}
	}//length是从0开始，所以+1
	if(cn_loc!=cn_length+1||diffa_loc!=diffa_length+1||diffb_loc!=diffb_length+1)
		{printf("error in find common neighbors(IN FUNTION FindNeighbors)\n");exit(1);}

	CN[0] = cn_length;
	DiffA[0] = diffa_length;
	DiffB[0] = diffb_length;

}


/* find CN */
void FindCN(int node_i, int node_j){
//node_1为节点1的独有邻居，node_j为节点2
	int num1 = NodeTable[node_i-1].NodeNum;
	int * A1 = NodeTable[node_i-1].Pointer;
	int num2 = NodeTable[node_j-1].NodeNum;
	int * A2 = NodeTable[node_j-1].Pointer;

	int p1_loc = 0;
	int p2_loc = 0;
	int * p1 = &A1[0];
	int * p2 = &A2[0];

	int cn_length = 0;
	int diffa_length = 0;
	int diffb_length = 0;

	int cn_loc = 1;
	int diffa_loc = 1;
	int diffb_loc = 1;

	while(p1_loc<num1 && p2_loc<num2){
		if(p1[p1_loc]<p2[p2_loc])
			p1_loc++;
		else if(p1[p1_loc]>p2[p2_loc])
			p2_loc++;
		else {
			CN[cn_loc] = p1[p1_loc];
			cn_length++;
			cn_loc++;
			p1_loc++;
			p2_loc++;
		}
	}

	if(cn_loc!=cn_length+1)
		{printf("error in find common neighbors(IN FUNTION FindCN)\n");exit(1);}
	CN[0] = cn_length;

}

/* sort each node's neighbors in ascending order 【以升序为每个节点的邻居排序】*/
void SortFun(int * Pointer, int Num){//Num为每个节点邻居节点数，指针从0开始
    //这里只是针对第i个节点
	int i = Num-1;//最大邻居
	int swap;
	while(i>0){
		int LastChangeIndex = 0;//最后改变位置
		for(int j=0;j<i;j++){
			if(Pointer[j]>Pointer[j+1]){//如果第j个指针所指向的节点id>第j+1个指向的节点id
				swap = Pointer[j+1];
				Pointer[j+1] = Pointer[j];
				Pointer[j] = swap;//交换指针所指向的节点
				LastChangeIndex = j+1;
			}
		}
		i = LastChangeIndex;
	}
}

int main(int argc, char * argv[]){
//  ./attractor Network_karate.txt Result_karate.txt 78IMPORTANT NOTES:
	// establish NodeTable
	printf("\nbegin establishing NodeTable...\n");
	EstablishNodeTable(argv[1],atoi(argv[3]));
	//用到了main函数中传的参数argv[]，改成手动输入比较好
	//argv[1]表示输入文件,argv[2]表示输出文件，atoi(argv[3])表示将字符串”边的数量“转化为整型数
    // argc是命令行总的参数个数  argv[]是argc个参数，其中第0个参数是程序的全名，以后的参数命令行后面跟的用户输入的参数
	printf("finish NodeTable\n");

	// establish EdgeTable
	printf("\nbegin establishing EdgeTable...\n");
	EstablishEdgeTable();
	printf("finish EdgeTable\n");


	// establish EdgesOfNodeTable
	printf("\nbegin establishing EdgesOfNodeTable...\n");
	EstablishEdgesOfNodeTable();
	printf("finish EdgesOfNodeTable\n");

	// interaction
	printf("\nbegin interacting...\n");
	Interaction(BETA, atoi(argv[3]));
	printf("finish interacting\n");

	// save result
	printf("\nbegin outputting...\n");
	FindClusters(argv[2]);
	printf("finish outputting\n");


	return 0;
}
