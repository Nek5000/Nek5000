/********************************types.h************************************
Author: Henry M. Tufo III

e-mail: hmt@@cs.brown.edu (no I didn't stutter on the <at> but lgrind ...)

sn-mail: Division of Applied Mathematics, Brown University,Providence, RI 02912

Last Modification: 5.2.95
*********************************types.h***********************************/



/********************************types.h************************************
NOTES ON USAGE: 

*********************************types.h***********************************/


/* data structures used to implement adjacency list */
typedef struct node  *adj_node_ptr;


typedef struct node{
   int elmt;
   double val;
   adj_node_ptr next;
 } adj_node;


adj_node_ptr init_adj_list(int);
void reset_adj_list(adj_node_ptr, int);
int insert_adj_list(adj_node_ptr, int, int, double);
void free_adj_list(adj_node_ptr, int);
void print_adj_list(adj_node_ptr, int);
int is_adj_list_symm(adj_node_ptr, int);
int *extract_od(adj_node_ptr adj_list, int N);




