
#ifndef _sparse_matrix_h
#define _sparse_matrix_h


void SMI_init();
void SMI_free_matrix(int);
int SMI_read_hb(char *,int);

int SMI_print_id(int);
adj_node_ptr SMI_get_adj(int mat);
matrix_id_ptr SMI_extract_matrix(int mat);


#endif
