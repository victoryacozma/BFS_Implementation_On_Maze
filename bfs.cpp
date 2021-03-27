#include <stdlib.h>
#include <string.h>
#include "bfs.h"

NodeT* newNode(int key) {
	NodeT* p = (NodeT*)malloc(sizeof(NodeT));
	p->cheie = key;
	p->numar_fii = 0;
	return p;
}

NodeT* transformare1(int* parinte, int nr)
{
	int k = 0;
	NodeT* v[100];
	for (int i = 0; i < nr; i++){
		v[i] = newNode(i);
	}

	for (int j = 0; j < nr; j++)
	{
		if (parinte[j] == -1)
		{
			k = j;
		}
		else
		{
			v[parinte[j]]->copii[v[parinte[j]]->numar_fii] = v[j];
			v[parinte[j]]->numar_fii++;
		}
	}
	return v[k];
}

void afisare(NodeT* r, Point* a, int k)
{
	if ( r!= NULL)
	{
		for (int i = 0; i < k; i++)
		{
			printf("\t");
		}
		printf("( %d, %d)\n", a[r->cheie].row, a[r->cheie].col);
		for (int i = 0; i < r->numar_fii; i++)
		{
			afisare(r->copii[i], a, k + 1);
		}
	}
}

int get_neighbors(const Grid *grid, Point p, Point neighb[])
{
    // TODO: fill the array neighb with the neighbors of the point p and return the number of neighbors
    // the point p will have at most 4 neighbors (up, down, left, right)
    // avoid the neighbors that are outside the grid limits or fall into a wall
    // note: the size of the array neighb is guaranteed to be at least 4
	int k = 0;
	int col = p.col;
	int row = p.row;
	if (grid->mat[row][col] == 0) {
		if (grid->mat[row][col - 1] == 0) {
			neighb[k].col = col - 1;
			neighb[k].row = row;
			k++;
		}
		if (grid->mat[row][col + 1] == 0) {
			neighb[k].col = col + 1;
			neighb[k].row = row;
			k++;

		}
		if (grid->mat[row - 1][col] == 0) {
			neighb[k].col = col;
			neighb[k].row = row - 1;
			k++;
		}
		if (grid->mat[row + 1][col] == 0) {
			neighb[k].col = col;
			neighb[k].row = row + 1;
			k++;
		}
		return k;
	}
	else return 0;
	
}

void grid_to_graph(const Grid *grid, Graph *graph)
{
    //we need to keep the nodes in a matrix, so we can easily refer to a position in the grid
    Node *nodes[MAX_ROWS][MAX_COLS];
    int i, j, k;
    Point neighb[4];

    //compute how many nodes we have and allocate each node
    graph->nrNodes = 0;
    for(i=0; i<grid->rows; ++i){
        for(j=0; j<grid->cols; ++j){
            if(grid->mat[i][j] == 0){
                nodes[i][j] = (Node*)malloc(sizeof(Node));
                memset(nodes[i][j], 0, sizeof(Node)); //initialize all fields with 0/NULL
                nodes[i][j]->position.row = i;
                nodes[i][j]->position.col = j;
                ++graph->nrNodes;
            }else{
                nodes[i][j] = NULL;
            }
        }
    }
    graph->v = (Node**)malloc(graph->nrNodes * sizeof(Node*));
    k = 0;
    for(i=0; i<grid->rows; ++i){
        for(j=0; j<grid->cols; ++j){
            if(nodes[i][j] != NULL){
                graph->v[k++] = nodes[i][j];
            }
        }
    }

    //compute the adjacency list for each node
    for(i=0; i<graph->nrNodes; ++i){
        graph->v[i]->adjSize = get_neighbors(grid, graph->v[i]->position, neighb);
        if(graph->v[i]->adjSize != 0){
            graph->v[i]->adj = (Node**)malloc(graph->v[i]->adjSize * sizeof(Node*));
            k = 0;
            for(j=0; j<graph->v[i]->adjSize; ++j){
                if( neighb[j].row >= 0 && neighb[j].row < grid->rows &&
                    neighb[j].col >= 0 && neighb[j].col < grid->cols &&
                    grid->mat[neighb[j].row][neighb[j].col] == 0){
                        graph->v[i]->adj[k++] = nodes[neighb[j].row][neighb[j].col];
                }
            }
            if(k < graph->v[i]->adjSize){
                //get_neighbors returned some invalid neighbors
                graph->v[i]->adjSize = k;
                graph->v[i]->adj = (Node**)realloc(graph->v[i]->adj, k * sizeof(Node*));
            }
        }
    }
}

void free_graph(Graph *graph)
{
    if(graph->v != NULL){
        for(int i=0; i<graph->nrNodes; ++i){
            if(graph->v[i] != NULL){
                if(graph->v[i]->adj != NULL){
                    free(graph->v[i]->adj);
                    graph->v[i]->adj = NULL;
                }
                graph->v[i]->adjSize = 0;
                free(graph->v[i]);
                graph->v[i] = NULL;
            }
        }
        free(graph->v);
        graph->v = NULL;
    }
    graph->nrNodes = 0;
}

Node* createNode()
{
	Node* p = (Node*)malloc(sizeof(Node));
	p->next = NULL;
	return p;
}

void enqueue(Node** head, Node** tail, Node* s)
{
	//TODO: insert the given key in the last position of the list given by head and tail;
	Node* p = s;
	if (*tail == NULL && *head == NULL)
	{
		(*head) = p;
		(*tail) = p;
	}
	else
	{
		(*tail)->next = p;
		*tail = p;
	}
}


Node* dequeue(Node** head, Node** tail)
{
	//TODO: delete first list element
	Node* r = NULL;
	if (*head != NULL)
	{
		//NodeT* r;
		r = (*head);
		(*head) = (*head)->next;

		if (*head == NULL)
		{
			*tail = NULL;
		}
	}
	return r;
}

bool goala(Node* first, Node* last)
{
	return ((first == NULL) && (last == NULL));

}

void bfs(Graph* graph, Node* s, Operation* op)
{
	// TOOD: implement the BFS algorithm on the graph, starting from the node s
	// at the end of the algorithm, every node reachable from s should have the color BLACK
	// for all the visited nodes, the minimum distance from s (dist) and the parent in the BFS tree should be set
	// for counting the number of operations, the optional op parameter is received
	// since op can be NULL (when we are calling the bfs for display purposes), you should check it before counting:
	// if(op != NULL) op->count();

	Node* first = NULL;
	Node* last = NULL;
	for (int i = 0; i < graph->nrNodes; i++) {
		graph->v[i]->color = COLOR_WHITE;
		graph->v[i]->dist = -1;
		graph->v[i]->parent = NULL;
	}
	s->color = COLOR_GRAY;
	s->dist = 0;
	s->parent = NULL;

	enqueue(&first, &last, s);

	while (!goala(first, last))
	{
		Node* p = dequeue(&first, &last);
		//printf("Nodul %d:  ", p->key);
		int col = p->position.col;
		int row = p->position.row;

		for (int i = 0; i < p->adjSize; i++) {
			if ((p->adj[i]->color) == COLOR_WHITE) {
				p->adj[i]->color = COLOR_GRAY;
				p->adj[i]->dist = p->dist + 1;
				p->adj[i]->parent = p;
				enqueue(&first, &last, p->adj[i]);
			}
		}
		p->color = COLOR_BLACK;
	}
}


void print_bfs_tree(Graph *graph)
{
    //first, we will represent the BFS tree as a parent array
    int n = 0; //the number of nodes
    int *p = NULL; //the parent array
    Point *repr = NULL; //the representation for each element in p

    //some of the nodes in graph->v may not have been reached by BFS
    //p and repr will contain only the reachable nodes
    int *transf = (int*)malloc(graph->nrNodes * sizeof(int));
    for(int i=0; i<graph->nrNodes; ++i){
        if(graph->v[i]->color == COLOR_BLACK){
            transf[i] = n;
            ++n;
        }else{
            transf[i] = -1;
        }
    }
    if(n == 0){
        //no BFS tree
        free(transf);
        return;
    }

    int err = 0;
    p = (int*)malloc(n * sizeof(int));
    repr = (Point*)malloc(n * sizeof(Node));
    for(int i=0; i<graph->nrNodes && !err; ++i){
        if(graph->v[i]->color == COLOR_BLACK){
            if(transf[i] < 0 || transf[i] >= n){
                err = 1;
            }else{
                repr[transf[i]] = graph->v[i]->position;
                if(graph->v[i]->parent == NULL){
                    p[transf[i]] = -1;
                }else{
                    err = 1;
                    for(int j=0; j<graph->nrNodes; ++j){
                        if(graph->v[i]->parent == graph->v[j]){
                            if(transf[j] >= 0 && transf[j] < n){
                                p[transf[i]] = transf[j];
                                err = 0;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    free(transf);
    transf = NULL;

    if(!err){
        // TODO: pretty print the BFS tree
        // the parrent array is p (p[k] is the parent for node k or -1 if k is the root)
        // when printing the node k, print repr[k] (it contains the row and column for that point)
        // you can adapt the code for transforming and printing multi-way trees from the previous labs
		NodeT* t = transformare1(p, n);
		afisare(t, repr, 0);
    }

    if(p != NULL){
        free(p);
        p = NULL;
    }
    if(repr != NULL){
        free(repr);
        repr = NULL;
    }
}

int shortest_path(Graph *graph, Node *start, Node *end, Node *path[])
{
    // TODO: compute the shortest path between the nodes start and end in the given graph
    // the nodes from the path, should be filled, in order, in the array path
    // the number of nodes filled in the path array should be returned
    // if end is not reachable from start, return -1
    // note: the size of the array path is guaranteed to be at least 1000
	bfs(graph, start);

	int m = end->dist;
	int k = 0;

	if (m == 0)
		return -1;
	while (start != end)
	{
		Node* n = end->parent;
		path[k] = n;
		k++;
		end = end->parent;
	}
	return m;
    
}

bool exista(int** mat, int a1, int a2) {

	
	if (mat[a1][a2] != 0)
			return true;
	return false;

}

void performance()
{
    int n, i;
    Profiler p("bfs");
	

    // vary the number of edges
    for(n=1000; n<=4500; n+=100){
        Operation op = p.createOperation("bfs-edges", n);
        Graph graph;
        graph.nrNodes = 100;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for(i=0; i<graph.nrNodes; ++i){
            graph.v[i] = (Node*)malloc(sizeof(Node));
			graph.v[i]->adjSize = graph.nrNodes;
			graph.v[i]->adj = (Node**)malloc(graph.v[i]->adjSize * sizeof(Node*));
            memset(graph.v[i], 0, sizeof(Node));
        }
        // TODO: generate n random edges
        // make sure the generated graph is connected
		int aux = 0;
		srand(time(NULL));
		int** mat = (int**)calloc(n * n, sizeof(int*));// cu arrayul acesta voi monitoriza ca adjunctii sa nu se repete
		
		//mai intai ne asiguram ca graful e conex
		//primul nod il vom lega cu toate celelalt
		
		for (int i = 0; i < graph.nrNodes; i++) {
			graph.v[0]->adj[i] = graph.v[i];
			printf("intra in for?");
			mat[graph.nrNodes][i] = 1;
			mat[i][graph.nrNodes] = 1;
			
		}
		//restul muchiilor
		for (int i = 0; i < n; i++) {
			int a1 = rand() % graph.nrNodes;
			int a2 = rand() % graph.nrNodes;

			printf("prima data se genereaza %d %d ", a1, a2);
				while ((mat[a1][a2] == 1) || (a1 == a2)) {
					printf("hello?");
					a1 = rand() % graph.nrNodes;
					a2 = rand() % graph.nrNodes;
					printf("apoi cel mai probabil aici crapa  cand cauta %d %d", a1, a2);
				}
				mat[a1][a2] = 1;
				mat[a2][a1] = 1;
				int aux1 = graph.v[a1]->adjSize;
				int aux2 = graph.v[a2]->adjSize;
				graph.v[a1]->adj[aux1] = graph.v[a2];
				graph.v[a2]->adj[aux2] = graph.v[a1];

				(graph.v[a1]->adjSize)++;
				(graph.v[a2]->adjSize)++;
		}
        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    // vary the number of vertices
    for(n=100; n<=200; n+=10){
        Operation op = p.createOperation("bfs-vertices", n);
        Graph graph;
        graph.nrNodes = n;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for(i=0; i<graph.nrNodes; ++i){
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }
        // TODO: generate 4500 random edges
        // make sure the generated graph is connected

		int aux = 0;
		srand(time(NULL));
		int** mat = (int**)calloc(n * n, sizeof(int*));// cu arrayul acesta voi monitoriza ca adjunctii sa nu se repete

		//mai intai ne asiguram ca graful e conex
		//primul nod il vom lega cu toate celelalte
		graph.v[0]->adjSize = graph.nrNodes;
		graph.v[0]->adj = (Node**)malloc(graph.v[i]->adjSize * sizeof(Node*));
		for (int i = 0; i < graph.nrNodes; i++) {
			graph.v[0]->adj[i] = graph.v[i];
			mat[graph.nrNodes][i] = 1;
			mat[i][graph.nrNodes] = 1;
			graph.v[i]->adjSize = 0;
		}
		//restul muchiilor
		for (int i = 0; i < 4500; i++) {
			int a1 = rand() % graph.nrNodes;
			int a2 = rand() % graph.nrNodes;

			printf("prima data se genereaza %d %d ", a1, a2);
			graph.v[i]->adj = (Node**)malloc(graph.v[i]->adjSize * sizeof(Node*));
			printf("se aloca matricea");
			while ((mat[a1][a2] == 1) || (a1 == a2)) {
				a1 = rand() % graph.nrNodes;
				a2 = rand() % graph.nrNodes;
				printf("apoi cel mai probabil aici crapa  cand cauta %d %d", a1, a2);
			}
			mat[a1][a2] = 1;
			mat[a2][a1] = 1;
			int aux1 = graph.v[a1]->adjSize;
			int aux2 = graph.v[a2]->adjSize;
			graph.v[a1]->adj[aux1] = graph.v[a2];
			graph.v[a2]->adj[aux2] = graph.v[a1];

			(graph.v[a1]->adjSize)++;
			(graph.v[a2]->adjSize)++;
		}

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    p.showReport();
}
