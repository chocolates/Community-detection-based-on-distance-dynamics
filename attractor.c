#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


// the Cohesion parameter $\lambda$ in our paper
#define BETA 0.6

// Network's node number
#define NetNodeNum 34

// suppose max neighbors number(the biggest degree) is not bigger than 1000
#define MaxNeighborNum 1000



// edge information
typedef struct Edge{
	int n_1_ID;
	int n_2_ID;
	int e_ID;
	double dist;

	int * CN_Loc;// CN_Loc points to a 2D matrix with 4 columns
	int CN_Loc_Length;// row number of CN_Loc, for search CN_Loc in Interaction
	
	int * N_1_Info;
	int N_1_Info_Length;// row number of N_1_Info, for search N_1_Info in Interaction

	int * N_1_Loc;
	int N_1_Loc_Length;// row number of N_1_Loc, for search N_1_Loc in Interaction

	int * N_2_Info;
	int N_2_Info_Length;// row number of N_2_Info, for search N_2_Info in Interaction

	int * N_2_Loc;
	int N_2_Loc_Length;// row number of N_2_Loc, for search N_2_Loc in Interaction

}EA;


// for EdgesOfNodeTable
typedef struct EdgeLoc{
	int Row;
	int Column;
}EdgeLoc;


// (NodeTable header)
typedef struct PointerToNeighbors{
	int NodeNum;// neighbors of the node
	int * Pointer;// points to the neighbors array
}NodeTableHead;


// (EdgeWithAttachmentTable header)
typedef struct PointerToEdges{
	int EdgeNum;
	EA * Pointer;
}EdgeTableHead;




// (EdgesOfNodeTable header)
typedef struct PointerToEdgeLoc{
	int EdgeNum;
	EdgeLoc * Pointer;
}EdgesOfNodeTableHead;


/* global variable */
NodeTableHead NodeTable[NetNodeNum]; // save the network, suppose node id begin at 1
EdgeTableHead EdgeWithAttachmentTable[NetNodeNum]; // save edges with its attachments
EdgesOfNodeTableHead EdgesOfNodeTable[NetNodeNum]; // save each node's edges(every edges that connected to it)
int CN[MaxNeighborNum];// most neighbors number is 66
int DiffA[MaxNeighborNum];
int DiffB[MaxNeighborNum];
int NetEdgeNum;
// other global variables: EdgesOfNodeTable, clusters

/* Establish NodeTable */
void EstablishNodeTable(char * inputfilename, int FileLine){
	FILE * file = fopen(inputfilename,"r");
	if(file==NULL)
		{printf("cannot open %s\n", inputfilename);exit(1);}
	else
		printf("opened %s\n successfully", inputfilename);
	int NeighborNumber[NetNodeNum];
	for(int i=0;i<NetNodeNum;i++)
		NeighborNumber[i] = 0;
	int node_1;
	int node_2;

	for(int i=0;i<FileLine;i++){
		fscanf(file,"%d %d",&node_1,&node_2);
		NeighborNumber[node_1-1]++;
		NeighborNumber[node_2-1]++;
	}

	// initiate NodeTable
	for(int i=0;i<NetNodeNum;i++){
		NodeTable[i].NodeNum = 0;
		NodeTable[i].Pointer = (int *)malloc(NeighborNumber[i]*sizeof(int));
		if(NodeTable[i].Pointer==NULL)
			{printf("Memory is not enough when initiating NodeTable\n");exit(1);}
	}
	fclose(file);
	file = fopen(inputfilename,"r");
	if(file==NULL)
		{printf("cannot open %s\n", inputfilename);exit(1);}
	else
		printf("have opened %s\n", inputfilename);

	for(int i=0;i<FileLine;i++){
		fscanf(file,"%d %d",&node_1,&node_2);
		NodeTable[node_1-1].Pointer[NodeTable[node_1-1].NodeNum++] = node_2;
		NodeTable[node_2-1].Pointer[NodeTable[node_2-1].NodeNum++] = node_1;
	}
	// check whether NodeTable[i].NodeNum equals to NeighborNumber[i] or not
	for(int i=0;i<NetNodeNum;i++)
		if(NodeTable[i].NodeNum!=NeighborNumber[i])
			{printf("NodeTable[%d] error\n", i+1);exit(1);}
	fclose(file);


	// sort each node's neighbors in ascending order
	void SortFun(int * , int );//declaration
	for(int i=0;i<NetNodeNum;i++)
		SortFun(NodeTable[i].Pointer, NodeTable[i].NodeNum);

}

/* Establish EdgeTable */
void EstablishEdgeTable(){
	int EdgeNumber[NetNodeNum];// save edge number of each node
	//compute node[loop_i]'s neighbors number
	for(int loop_i=0;loop_i<NetNodeNum;loop_i++){
		EdgeNumber[loop_i] = 0;
		for(int loop_j=0; loop_j<NodeTable[loop_i].NodeNum; loop_j++)
			if(NodeTable[loop_i].Pointer[loop_j]>loop_i+1)// save each edge only once
				EdgeNumber[loop_i]++;
	}

	// compute n_1_ID, n_2_ID, e_ID
	int edgeid = 0;
	for(int loop_i=0;loop_i<NetNodeNum;loop_i++){
		EdgeWithAttachmentTable[loop_i].EdgeNum = EdgeNumber[loop_i];
		EdgeWithAttachmentTable[loop_i].Pointer = (EA *)malloc(EdgeNumber[loop_i]*sizeof(EA));
		if(EdgeWithAttachmentTable[loop_i].Pointer==NULL)
			{printf("memory is not enough with node%d's edges\n", loop_i+1);exit(1);}

		int edge_loc_tmp=0;
		for(int loop_j=0; loop_j<NodeTable[loop_i].NodeNum; loop_j++){
			if(NodeTable[loop_i].Pointer[loop_j]>loop_i+1){
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].n_1_ID = loop_i+1;
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].n_2_ID = NodeTable[loop_i].Pointer[loop_j];
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].e_ID = edgeid;
				edgeid++;
				NetEdgeNum = edgeid;
				edge_loc_tmp++;
			}
		}
		// check whether edge_loc_tmp equals to EdgeWithAttachmentTable[loop_i].EdgeNum or not
		if(edge_loc_tmp!=EdgeWithAttachmentTable[loop_i].EdgeNum)
			{printf("EdgeWithAttachmentTable[%d] error(line 175)\n", loop_i+1); exit(1);}
	}

	/* compute dist, CN_Loc, N_1_Info, N_1_Loc, N_2_Info, N_2_Loc */
	for(int i=0;i<NetNodeNum;i++){
		for(int j=0;j<EdgeWithAttachmentTable[i].EdgeNum;j++){
			int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
			int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
			void FindNeighbors(int, int);// function declaration
			FindNeighbors(node_1,node_2);
			
			// compute distance
			EdgeWithAttachmentTable[i].Pointer[j].dist = 1.0 - (double)(CN[0]+2)/(double)(CN[0]+DiffA[0]+DiffB[0]+2);
			
			/* establish CN_Loc */
			EdgeWithAttachmentTable[i].Pointer[j].CN_Loc = (int *)malloc(4*CN[0]*sizeof(int));
			if(EdgeWithAttachmentTable[i].Pointer[j].CN_Loc==NULL)
				{printf("EdgeWithAttachmentTable[%d].Pointer is NULL\n", i);exit(1);}
			EdgeWithAttachmentTable[i].Pointer[j].CN_Loc_Length = CN[0];
			for(int s=1;s<=CN[0];s++){
				int NodeCN = CN[s];
				int Loc_Node_1_R = -1;
				int Loc_Node_1_C = -1;
				int Loc_Node_2_R = -1;
				int Loc_Node_2_C = -1;
				int NodeMin = (node_1<NodeCN)?node_1:NodeCN;
				int NodeMax = (node_1>NodeCN)?node_1:NodeCN;
				Loc_Node_1_R = NodeMin - 1;
				for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
					if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
						{	Loc_Node_1_C = loop;	break;}
				if(Loc_Node_1_C==-1)
					{printf("line 207 error\n"); exit(1);}
				NodeMin = (node_2<NodeCN)?node_2:NodeCN;
				NodeMax = (node_2>NodeCN)?node_2:NodeCN;
				Loc_Node_2_R = NodeMin - 1;
				for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
					if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
						{	Loc_Node_2_C = loop;	break;}
				if(Loc_Node_2_C==-1)
					{printf("line 261 error\n"); exit(1);}
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4*(s-1)+0] = Loc_Node_1_R;
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4*(s-1)+1] = Loc_Node_1_C;
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4*(s-1)+2] = Loc_Node_2_R;
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4*(s-1)+3] = Loc_Node_2_C;
			}

//			printf("Node %d, Node %d complete CN_Loc\n", node_1, node_2);


			/* establish N_1_Loc, N_1_Info */

			EdgeWithAttachmentTable[i].Pointer[j].N_1_Info = (int *)malloc(2*2*DiffA[0]*sizeof(int));// each node have four elements: node id, the number common neighbors of N_1 and node_2, row of <N_1, node_1>, column of <N_1, node_1>
			EdgeWithAttachmentTable[i].Pointer[j].N_1_Info_Length = DiffA[0];

			if(EdgeWithAttachmentTable[i].Pointer[j].N_1_Info==NULL)
				{printf("memory not enough With N_1_Info\n"); exit(1);}
			int Edgenum_tmp = 0;
			for(int s=1;s<=DiffA[0];s++){
				int Node_N_1 = DiffA[s];
				void FindCN(int, int);// function declaration
				FindCN(Node_N_1, node_2);// CN changes
				Edgenum_tmp  = Edgenum_tmp + CN[0];
			}	
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
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+0] = Node_N_1;// this node's id
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+1] = CN[0];// number or common neighbors of this node and node_2

				int NodeMin;
				int NodeMax;

				// N_1_Info[4*(s-1)+2] is the location of row of <N_1, node_1>, N_1_Info[4*(s-1)+3] is the column of <N_1, node_1>
				NodeMin = (node_1<Node_N_1)?node_1:Node_N_1;
				NodeMax = (node_1>Node_N_1)?node_1:Node_N_1;
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+2] = NodeMin-1;//row
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+3] = -1;
				for(int loop=0;loop<EdgeWithAttachmentTable[NodeMin-1].EdgeNum;loop++)
					if(EdgeWithAttachmentTable[NodeMin-1].Pointer[loop].n_2_ID == NodeMax)
						{	EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4*(s-1)+3] = loop;	break;}
				if(EdgeWithAttachmentTable[NodeMin-1].EdgeNum==-1)
					{printf("line 270 error\n");exit(1);}	
					

				for(int ss=1;ss<=CN[0];ss++){
					int NodeCN = CN[ss];
					int Loc_Node_N_1_R = -1;
					int Loc_Node_N_1_C = -1;
					int Loc_Node_2_R = -1;
					int Loc_Node_2_C = -1;

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

			/* establish N_2_Loc, N_2_Info */
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
		EdgesOfNodeTable[i].EdgeNum = 0;
		EdgesOfNodeTable[i].Pointer = (EdgeLoc *)malloc(NodeTable[i].NodeNum*sizeof(EdgeLoc));
		if(EdgesOfNodeTable[i].Pointer==NULL)
			{printf("memory not enough with EdgesOfNodeTable(line 378)\n");exit(1);}
	}
	for(int i=0;i<NetNodeNum;i++)
		for(int j=0;j<EdgeWithAttachmentTable[i].EdgeNum;j++){
			int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
			int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
			// a little check
			if(node_1!=i+1)
				{printf("line 386 error\n"); exit(1);}

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

/* Interaction */
void Interaction(double beta, int FileLine){
	// beta controls clusters number
	if(NetEdgeNum!=FileLine){
		printf("NetEdgeNum error...(line 413)\n");
		exit(1);
	}
	double D[NetEdgeNum];

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
	int Terminate = 1;
	int Loop = 0;
	while(Terminate){
		Loop++;
		printf("Loop: %d\n", Loop);
		for(int s_1=0;s_1<NetNodeNum;s_1++)
			for(int s_2=0;s_2<EdgeWithAttachmentTable[s_1].EdgeNum;s_2++)
				if(EdgeWithAttachmentTable[s_1].Pointer[s_2].dist>0 && EdgeWithAttachmentTable[s_1].Pointer[s_2].dist<1){
					EA ThisEA = EdgeWithAttachmentTable[s_1].Pointer[s_2];
					int ThisNode_1 = ThisEA.n_1_ID;
					int ThisNode_2 = ThisEA.n_2_ID;
                    int ThisEdgeId = ThisEA.e_ID;

					double CI = 0.0;// common neighbors' influence
					double N_1_I = 0.0;// node_1's neighbors' influence
					double N_2_I = 0.0;// node_2's neighbors' influence

					/* common neighbors' influence */
					for(int s_3=0;s_3<ThisEA.CN_Loc_Length;s_3++){
						double Distance_CI_1;
						double Distance_CI_2;
						int * CNLocTmp = ThisEA.CN_Loc;

						Distance_CI_1 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+0]].Pointer[CNLocTmp[s_3*4+1]].dist;
						Distance_CI_2 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+2]].Pointer[CNLocTmp[s_3*4+3]].dist;

						int CNTmp_1 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+0]].Pointer[CNLocTmp[s_3*4+1]].n_1_ID;
						int CNTmp_2 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+0]].Pointer[CNLocTmp[s_3*4+1]].n_2_ID;
						int CNTmp_3 = (CNTmp_1==ThisNode_1)?CNTmp_2:CNTmp_1;
						CI = CI + ( (double)(1-Distance_CI_2) / (double)NodeTable[ThisNode_1-1].NodeNum ) *  sin(1-Distance_CI_1);
						// printf("CI influence1: %f\n", ( (double)(1-Distance_CI_2) / (double)NodeTable[ThisNode_1-1].NodeNum ) *  sin(1-Distance_CI_1));
						
						CNTmp_1 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+2]].Pointer[CNLocTmp[s_3*4+3]].n_1_ID;
						CNTmp_2 = EdgeWithAttachmentTable[CNLocTmp[s_3*4+2]].Pointer[CNLocTmp[s_3*4+3]].n_2_ID;
						CNTmp_3 = (CNTmp_1==ThisNode_2)?CNTmp_2:CNTmp_1;
						CI = CI + ( (double)(1-Distance_CI_1) / (double)NodeTable[ThisNode_2-1].NodeNum ) *  sin(1-Distance_CI_2);
						// printf("CI influence2: %f\n", ( (double)(1-Distance_CI_1) / (double)NodeTable[ThisNode_2-1].NodeNum ) *  sin(1-Distance_CI_2));
					}

					/* node_1's neighbors' influence */
					int s_3 = 0;
					int s_4 = 0;
					int s_4_count;
					while(s_3<ThisEA.N_1_Info_Length){
						int Ngh_N_1 = ThisEA.N_1_Info[s_3*4+0];// node_1's neighbor's id
						int CN_Num = ThisEA.N_1_Info[s_3*4+1];// common neighbors' number of node_1's neighbor and node_2
						double Distance_N_N_1 = EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3*4+2]].Pointer[ThisEA.N_1_Info[s_3*4+3]].dist;

						s_4_count = s_4;
						double lambda_numerator = 0;
						for(s_4=s_4_count;s_4<s_4_count+CN_Num;s_4++){
							int loc_1;
							int loc_2;

							loc_1 = ThisEA.N_1_Loc[s_4*4+0];
							loc_2 = ThisEA.N_1_Loc[s_4*4+1];
							lambda_numerator = lambda_numerator + (1-EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
							loc_1 = ThisEA.N_1_Loc[s_4*4+2];
							loc_2 = ThisEA.N_1_Loc[s_4*4+3];
							lambda_numerator = lambda_numerator + (1-EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
						}

						double lambda_denominator = 0;
						for(int s_5=0;s_5<EdgesOfNodeTable[Ngh_N_1-1].EdgeNum;s_5++){
							int loc_r = EdgesOfNodeTable[Ngh_N_1-1].Pointer[s_5].Row;
							int loc_c = EdgesOfNodeTable[Ngh_N_1-1].Pointer[s_5].Column;
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
		if(sum_1==sum_2 || Loop>1000)//end condition
			Terminate = 0;

		// renew edge's distance
		EdgeCounter = 0;
		for(int s_1=0;s_1<NetNodeNum;s_1++)
			for(int s_2=0;s_2<EdgeWithAttachmentTable[s_1].EdgeNum;s_2++){
				EdgeWithAttachmentTable[s_1].Pointer[s_2].dist = D[EdgeCounter];
				EdgeCounter++;
			}

	}

}


/* find clusters, by BFS */
void FindClusters(char * outputfile){
	int CountNeighborNumber[NetNodeNum];
	// establish CountNeighborNumber
	for(int i=0;i<NetNodeNum;i++)
		CountNeighborNumber[i] = 0;
	for(int i=0;i<NetNodeNum;i++)
		for(int j=0;j<EdgeWithAttachmentTable[i].EdgeNum;j++)
			if(EdgeWithAttachmentTable[i].Pointer[j].dist<1){
				int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
				int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
				CountNeighborNumber[node_1-1]++;
				CountNeighborNumber[node_2-1]++;
			}
	// establish new NodeTable
	NodeTableHead NodeTable_New[NetNodeNum];
	for(int i=0;i<NetNodeNum;i++){
		NodeTable_New[i].NodeNum = 0;
		NodeTable_New[i].Pointer = (int *)malloc(CountNeighborNumber[i]*sizeof(int));
		if(NodeTable_New[i].Pointer==NULL){
			printf("memory is not enough with NodeTable_New\n");
			exit(1);
		}
	}
	for(int i=0;i<NetNodeNum;i++)
		for(int j=0;j<EdgeWithAttachmentTable[i].EdgeNum;j++)
			if(EdgeWithAttachmentTable[i].Pointer[j].dist<1){
				int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
				int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
				NodeTable_New[node_1-1].Pointer[NodeTable_New[node_1-1].NodeNum] = node_2;
				NodeTable_New[node_2-1].Pointer[NodeTable_New[node_2-1].NodeNum] = node_1;
				NodeTable_New[node_1-1].NodeNum++;
				NodeTable_New[node_2-1].NodeNum++;
			}

	//check if NodeTable is ok
	for(int i=0;i<NetNodeNum;i++)
		if(NodeTable_New[i].NodeNum!=CountNeighborNumber[i])
			printf("NodeTable_New Error! In Node: %d\n", i+1);

	/* find clusters */
	int clusters[NetNodeNum+1];// clusters[i] is node_i's id, clusters[0] is unused
	for(int k=0;k<=NetNodeNum;k++)
		clusters[k] = -1;
	int Q[NetNodeNum];// the queue
	int ClusterID = 1;
	int Terminate = 1;

	while(Terminate){
		Terminate = 0;
		int id;
		for(id=1;id<=NetNodeNum;id++)
			if(clusters[id]==-1){
				Terminate = 1;
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

	// save clusters in 2D array
	FILE * fout = fopen(outputfile,"w");
	if(fout==NULL)
		{	printf("opening outputfile fails\n"); exit(0);}
	for(int i=1;i<=NetNodeNum;i++)
		fprintf(fout, "%d %d\n", i, clusters[i]);
	fclose(fout);
	


}

/* compute cn, diffA, diffB */
void FindNeighbors(int node_i, int node_j){

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
			CN[cn_loc] = A1[p1_loc];
			cn_length++;
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
	}
	if(cn_loc!=cn_length+1||diffa_loc!=diffa_length+1||diffb_loc!=diffb_length+1)
		{printf("error in find common neighbors(IN FUNTION FindNeighbors)\n");exit(1);}

	CN[0] = cn_length;
	DiffA[0] = diffa_length;
	DiffB[0] = diffb_length;

}


/* find CN */
void FindCN(int node_i, int node_j){

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

/* sort each node's neighbors in ascending order */
void SortFun(int * Pointer, int Num){
	int i = Num-1;
	int swap;
	while(i>0){
		int LastChangeIndex = 0;
		for(int j=0;j<i;j++){
			if(Pointer[j]>Pointer[j+1]){
				swap = Pointer[j+1];
				Pointer[j+1] = Pointer[j];
				Pointer[j] = swap;
				LastChangeIndex = j+1;
			}
		}
		i = LastChangeIndex;
	}
}

int main(int argc, char * argv[]){

	// establish NodeTable
	printf("\nbegin establishing NodeTable...\n");
	EstablishNodeTable(argv[1],atoi(argv[3]));
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





