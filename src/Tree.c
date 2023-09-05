#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Utils.h"
#include "Tree.h"


/*print ident, time and comment of node n time_name*/	
void fprintIdentTimeComment(FILE *f, int n, TypeTree *tree, TypeDisplayName display) {
	switch(display) {
		case display_none:
		case display_time_none:
			break;
		case display_name:
		case display_time_name:
			if(tree->name != NULL && tree->name[n] != NULL)
				fprintf(f, "'%s'", tree->name[n]);
			break;
		case display_index:
		case display_time_index:
			fprintf(f, "'%d'", n);
			break;
		case display_both:
		case display_time_both:
		default:
			if(tree->name != NULL && tree->name[n] != NULL)
				fprintf(f, "'%s-", tree->name[n]);
			else
				fprintf(f, "'");
			fprintf(f, "%d'", n);
			break;
	}
	if(display>=display_time_none &&  tree->time != NULL)
		fprintf(f, ":%lf", tree->time[n]);
	if(tree->comment != NULL && tree->comment[n] != NULL)
		fprintf(f, "[%s]", tree->comment[n]);
}

/*print tree in newick format*/
void fprintTreeNewick(FILE *f, TypeTree *tree) {
    if(tree->size<=0)
        return;
    if(tree->node[tree->root].child >= 0) {
        int tmp = tree->node[tree->root].child;
        fprintf(f, "(");
        fprintNodeNewick(f, tmp, tree);
        for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
            fprintf(f, ", ");
            fprintNodeNewick(f, tmp, tree);
        }
        fprintf(f, ")");
    }
    fprintIdentTimeComment(f, tree->root, tree, display_time_name);
    fprintf(f, ";\n");
}

/*print node in newick format*/
void fprintNodeNewick(FILE *f, int n, TypeTree *tree) {
    if(tree->node[n].child >= 0) {
        int tmp = tree->node[n].child;
        fprintf(f, "(");
        fprintNodeNewick(f, tmp, tree);
        for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
            fprintf(f, ", ");
            fprintNodeNewick(f, tmp, tree);
        }
        fprintf(f, ")");
    }
    fprintIdentTimeComment(f, n, tree, display_time_name);
}


/*fully duplicate "tree"*/
TypeTree *cpyTree(TypeTree *tree) {
    int n;
    TypeTree *res;

    res = (TypeTree*) malloc(sizeof(TypeTree));
    res->sizeBuf = tree->sizeBuf;
    res->size = tree->size;
    res->node = (TypeNode*) malloc(res->sizeBuf*sizeof(TypeNode));
    for(n=0; n<tree->sizeBuf; n++)
        res->node[n] = tree->node[n];
    if(tree->time) {
        res->time = (double*) malloc(res->sizeBuf*sizeof(double));
        for(n=0; n<tree->sizeBuf; n++)
            res->time[n] = tree->time[n];
    } else
        res->time = NULL;
    if(tree->parent) {
        res->parent = (int*) malloc(res->sizeBuf*sizeof(int));
        for(n=0; n<tree->sizeBuf; n++)
            res->parent[n] = tree->parent[n];
    } else
        res->parent = NULL;
    if(tree->name) {
        res->name = (char**) malloc(res->sizeBuf*sizeof(char*));
        for(n=0; n<tree->sizeBuf; n++)
            if(tree->name[n] != NULL)
                res->name[n] = strdpl(tree->name[n]);
            else
                res->name[n] = NULL;
    } else
        res->name = NULL;
    if(tree->comment) {
        res->comment = (char**) malloc(res->sizeBuf*sizeof(char*));
        for(n=0; n<tree->sizeBuf; n++)
            if(tree->comment[n] != NULL)
                res->comment[n] = strdpl(tree->comment[n]);
            else
                res->comment[n] = NULL;
    } else
        res->comment = NULL;
    res->info = NULL;
    res->root = tree->root;
    return res;
}


void fillParent(int node, TypeTree *tree, int *parent) {
	int c;
	for(c=tree->node[node].child; c!=NOSUCH; c=tree->node[c].sibling) {
		fillParent(c, tree, parent);
		parent[c] = node;
	}
}

/*return the table of parents*/
int *getParent(TypeTree *tree) {
	int n, *parent;
	if(tree->size == 0)
		return NULL;
	parent = (int*) malloc(tree->sizeBuf*sizeof(int));
	for(n=0; n<tree->sizeBuf; n++)
		parent[n] = NOSUCH;
	fillParent(tree->root, tree, parent);
	return parent;
}


/*desallocate tree*/
void freeTree(TypeTree *tree) {
    if(tree == NULL)
        return;
    if(tree->node != NULL)
        free((void*)tree->node);
    if(tree->time != NULL)
        free((void*)tree->time);
    if(tree->parent != NULL)
        free((void*)tree->parent);
    if(tree->name != NULL) {
        int i;
        for(i=0; i<tree->sizeBuf; i++)
            if(tree->name[i] != NULL)
                free((void*)tree->name[i]);
        free((void*)tree->name);
    }
    if(tree->comment != NULL) {
        int i;
        for(i=0; i<tree->sizeBuf; i++)
            if(tree->comment[i] != NULL)
                free((void*)tree->comment[i]);
        free((void*)tree->comment);
    }
    free((void*)tree);
}

void fillTipsRec(int n, TypeTree *tree, int *tips, int *ntips) {
	if(tree->node[n].child == NOSUCH)
		tips[(*ntips)++] = n;
	else {
		int c;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling)
			fillTipsRec(c, tree, tips, ntips);
	}
}
	

int fillTips(int n, TypeTree *tree, int *tips) {
	int ntips = 0;
	fillTipsRec(n, tree, tips, &ntips);
	return ntips;
}


TypeTree *rerootTreeREML(int node, TypeTree *tree) {
	TypeTree *res;
	int *parent, n, *path, size, i, *toSetBis;
	if(tree->size<3)
		return NULL;
	parent = getParent(tree);
	path = (int*) malloc((tree->size+1)*sizeof(int));
	size=0;
	for(n=node; n!=NOSUCH; n=parent[n])
		path[size++] = n;
	path[size] = NOSUCH;
	res = cpyTree(tree);

	for(i=size-1; i>0; i--) {
		for(toSetBis=&(res->node[path[i]].child); *toSetBis!=NOSUCH && *toSetBis!=path[i-1];  toSetBis = &(res->node[*toSetBis].sibling))
			;
		if(*toSetBis!=NOSUCH)
			*toSetBis = res->node[*toSetBis].sibling;
		res->node[path[i]].sibling = res->node[path[i]].child;
		res->node[path[i]].child = parent[path[i]];
	}
	res->node[node].sibling = NOSUCH;
	res->node[node].child = parent[node];
	res->root = node;
	for(toSetBis=&(res->node[path[size-2]].child); *toSetBis!=NOSUCH && *toSetBis!=path[size-1];  toSetBis = &(res->node[*toSetBis].sibling))
		;
	if(*toSetBis!=NOSUCH)
		*toSetBis = res->node[*toSetBis].sibling;
	if(res->node[res->node[path[size-2]].child].sibling == NOSUCH) {
		if(size == 2) {
			res->time[path[size-2]] += res->time[res->node[path[size-2]].child];
			res->node[path[size-2]].child = res->node[res->node[path[size-2]].child].child;
		} else {
			for(toSetBis=&(res->node[path[size-3]].child); *toSetBis!=NOSUCH && *toSetBis!=path[size-2];  toSetBis = &(res->node[*toSetBis].sibling))
				;
			if(*toSetBis!=NOSUCH)
				*toSetBis = res->node[*toSetBis].sibling;
			res->time[res->node[path[size-2]].child] += res->time[path[size-2]];
			res->node[res->node[path[size-2]].child].sibling = res->node[path[size-3]].child;
			res->node[path[size-3]].child = res->node[path[size-2]].child;
		}
	}
	free((void*)parent);
	free((void*)path);
	res->size -= 2;
	return res;
}

TypeTree *rerootTreeStem(int node, TypeTree *tree) {
	TypeTree *res;
	int n, *parent, *path, size, i;
	if(tree->node[node].child != NOSUCH)
		return NULL;
	parent = getParent(tree);
	parent[tree->root] = NOSUCH;
	path = (int*) malloc(tree->size*sizeof(int));
	size=0;
	for(n=node; n!=NOSUCH; n=parent[n])
		path[size++] = n;
	res = cpyTree(tree);
	if(size<2)
		return res;
	int *toSetBis;		
	for(i=size-1; i>0; i--) {
		for(toSetBis=&(res->node[path[i]].child); *toSetBis!=NOSUCH && *toSetBis!=path[i-1];  toSetBis = &(res->node[*toSetBis].sibling))
			;
		if(*toSetBis!=NOSUCH)
			*toSetBis = res->node[*toSetBis].sibling;
		res->node[path[i]].sibling = res->node[path[i]].child;
		res->node[path[i]].child = parent[path[i]];
	}
	res->node[node].sibling = NOSUCH;
	res->node[node].child = parent[node];
	free((void*)parent);
	free((void*)path);
	res->root = node;
	return res;
}

int findSide(int node, TypeTree *tree) {
	int n, *parent;
	parent = getParent(tree);
	if(parent[node] == NOSUCH)
		return node;
	for(n=node; parent[parent[n]] != NOSUCH; n = parent[n])
		;
	free((void*)parent);
	return n;
}
