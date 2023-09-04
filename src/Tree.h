#ifndef TreeF
#define TreeF
#include <stdio.h>

#define NOSUCH -1
#define NO_TIME -DBL_MAX


typedef enum DISPLAY_NAME {
	display_none=0,
	display_name,
	display_index,
	display_both,
	display_time_none,
	display_time_name,
	display_time_index,
	display_time_both
} TypeDisplayName;


typedef struct NODE {
    int child, sibling;
} TypeNode;

typedef struct TREE {
    TypeNode *node;
    int root, size, sizeBuf, *parent;
    double *time;
    char **name, **comment;
    void *info;
} TypeTree;


#ifdef __cplusplus
extern "C" {
namespace Tree {
#endif

int fillTips(int n, TypeTree *tree, int *tips);
int *getParent(TypeTree *tree);
/*fully duplicate "tree"*/
TypeTree *cpyTree(TypeTree *tree);
/*desallocate tree*/
void freeTree(TypeTree *tree);
void fprintIdentTimeComment(FILE *f, int n, TypeTree *tree, TypeDisplayName display);
/*print tree in newick format*/
void fprintTreeNewick(FILE *f, TypeTree *tree);
/*print node in newick format*/
void fprintNodeNewick(FILE *f, int n, TypeTree *tree);
TypeTree *rerootTreeREML(int node, TypeTree *tree);
TypeTree *rerootTreeStem(int node, TypeTree *tree);
int findSide(int node, TypeTree *tree);

#ifdef __cplusplus
}
}
#endif

#endif
