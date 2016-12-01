# Attractor -- Community detection based on distance dynamics

We implement our proposed community detection algorithm (Attractor) in 'attractor.c'. 

## How to run the code:
1. change the global variable value(it is the cohesion parameter $\lambda$ in our paper) in line 8;
2. change the global variable value(it is the network's node number) in line 17;
3. Compile the codes;
	For example: On Mac: gcc attractor.c -o attractor;
			 On linux: gcc -std=c99 attractor.c -lm -o attractor
4. execute it with three input and one output: first input is the network file(its format is explained below), second input is the result file(its format is explained below), last input is the edge number;
	For example: In bash: ./attractor Network_karate.txt Result_karate.txt 78


## IMPORTANT NOTES:
- Format of the input network file(see ./Network_example.txt):
	1) each line, represents an undirected edge, contains two numbers: source node and the target node
	2) each edge is just saved only once;
	3) the numbers of node begin at 1 and they must be continuous, that is, if the network has 34 nodes, the 
	smallest node number is 1 and the biggest node number is 34

- Format of the output result(see ./Result_example.txt): each line contains two number: the node number and its corresponding cluster

- This source code supposes that the degree of the biggest degree node is smaller than 1000. If it is not this situation, the codes in line 12 should be modified according the real situation.

## UPDATE
### 2016-11-30
- Correct _SortFun_ function which may lead to mistakes in finding CN and EN. For experiments, results of Zarachy, Football, Polbooks, Amazon, Friendship and synthetic networks are not influenced. However, hepth (Collaboration) and RoadPa are influenced. HEPTH: correct modularity is 0.337 and new NCUT is 214.226 with 785 clusters (vs 0.579 and 1179 with 1384 clusters); ROADPA: correct modularity is 0.865 and new NCUT is 23136 with 56967 clusters (vs 0.856 and 25055 with 59919 clusters).