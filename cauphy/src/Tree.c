#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Tree.h"


#define INFTY 9.9E99;
#define BASIC_TMP_SIZE 1000
#define BASIC_INC_BUFFER 20


typedef struct TMP_REORDER {
    char *name;
    int index;
} TypeTmpReorder;

/*turn branch lengthes to absolute time rec*/
static void bltoabsTimeRec(int n, double offset, TypeTree *tree);
/*turn branch lengthes to absolute time rec*/
static void abstoblTimeRec(int n, double offset, TypeTree *tree);
static int compareTmpReorder(const void *a, const void *b);
static void reorderTreeRec(char **tmp, char **name, int n, TypeTree *tree);
static char *nameInternalNodesRec(int n, char **name, TypeTree *tree);
static int iterateBinary(int n, TypeTree *resT, TypeTree *tree);
static void printNodeDebug(FILE *f, int s, int depth, TypeTree *tree, char **name);
static void setAbsentRec(int n, int *index, TypeTree *tree);
static void fillTipsRec(int n, TypeTree *tree, int *tips, int *ntips);


void completeName(TypeTree *tree) {
	int n;
	if(tree->name == NULL) {
		tree->name = (char**) malloc(tree->sizeBuf*sizeof(char*));
		for(n=0; n<tree->sizeBuf; n++) {
			char tmp[100];
			sprintf(tmp, "%d", n);
			tree->name[n] = strdpl(tmp);
		}
	} else
		for(n=0; n<tree->size; n++) {
			if(tree->name[n] == NULL) {
				char tmp[100];
				sprintf(tmp, "%d", n);
				tree->name[n] = strdpl(tmp);
			}
		}
}

double getMaximumLeafTime(TypeTree *tree) {
    int n;
    for(n=0; n<tree->size && (tree->node[n].child != NOSUCH || tree->time[n] == NO_TIME); n++);
    if(n>=tree->size)
        return NO_TIME;
    else {
        double max = tree->time[n];
        n++;
        for(; n<tree->size; n++)
            if(tree->node[n].child == NOSUCH && tree->time[n] != NO_TIME && tree->time[n]>max)
                max = tree->time[n];
        return max;
    }
}


void setAbsentRec(int n, int *index, TypeTree *tree) {
    if(n == NOSUCH)
        return;
    index[n] = NOSUCH;
    int c;
    for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
        setAbsentRec(c, index, tree);
}

/*remove the subtree pending from n*/
int *removeSubtreeReturnIndex(int n, TypeTree *tree) {
    int *index, i, c, nchild;

    if(tree->parent == NULL)
        tree->parent = getParent(tree);

    index = (int*) malloc(tree->size*sizeof(int));
    for(i=0; i<tree->size; i++)
        index[i] = 1;
    setAbsentRec(n, index, tree);
    if(tree->parent[n] != NOSUCH) {
        for(c=tree->node[tree->parent[n]].child,nchild=0; c!=NOSUCH; c=tree->node[c].sibling,nchild++)
            ;
        if(nchild>2) {
            int *toSet;
            for(toSet=&(tree->node[tree->parent[n]].child); *toSet!=NOSUCH && *toSet!=n;  toSet = &(tree->node[*toSet].sibling))
                ;
            if(*toSet != n) {
                error("Issue A while removing  node %d (parent %d)\n", n, tree->parent[n]);
            } else
                *toSet = tree->node[n].sibling;
        } else { //remove parent[n]
            index[tree->parent[n]] = NOSUCH;
            for(c=tree->node[tree->parent[n]].child; c!=NOSUCH && c==n; c=tree->node[c].sibling)
                ;
            if(c==NOSUCH) {
                error("Issue B while removing  node %d (parent %d) - node with only one child\n", n, tree->parent[n]);
            }
            if(tree->parent[tree->parent[n]] == NOSUCH) {
                tree->root = c;
            } else {
                int *toSet;
                for(toSet=&(tree->node[tree->parent[tree->parent[n]]].child); *toSet!=NOSUCH && *toSet!=tree->parent[n];  toSet = &(tree->node[*toSet].sibling))
                    ;
                if(*toSet != tree->parent[n]) {
                    error("Issue C while removing  node %d (parent %d)\n", n, tree->parent[n]);
                } else {
                    tree->node[c].sibling = tree->node[tree->parent[n]].sibling;
                    *toSet = c;
                     tree->parent[c] = tree->parent[tree->parent[n]];
                }
            }
        }
    } else { //suppress root - nothing to do
    }
    int ind=0;
    for(i=0; i<tree->size; i++)
        if(index[i] != NOSUCH)
            index[i] = ind++;
        else {
             if(tree->name != NULL && tree->name[i] != NULL) {
                free((void*)tree->name[i]);
                tree->name[i] = NULL;
            }
            if(tree->comment != NULL && tree->comment[i] != NULL) {
                free((void*)tree->comment[i]);
                tree->comment[i] = NULL;
            }
        }
    for(i=0; i<tree->size; i++) {
        if(index[i] != NOSUCH) {
            if(tree->node[i].child != NOSUCH)
                tree->node[index[i]].child = index[tree->node[i].child];
            else
                tree->node[index[i]].child = NOSUCH;
            if(tree->node[i].sibling != NOSUCH)
                tree->node[index[i]].sibling = index[tree->node[i].sibling];
            else
                tree->node[index[i]].sibling = NOSUCH;
            if(tree->time != NULL)
                tree->time[index[i]] = tree->time[i];
            if(tree->name != NULL)
                tree->name[index[i]] = tree->name[i];
            if(tree->comment != NULL)
                tree->comment[index[i]] = tree->comment[i];
            if(tree->parent[i] != NOSUCH)
                tree->parent[index[i]] = index[tree->parent[i]];
            else
                tree->parent[index[i]] = NOSUCH;
        }
    }
    for(i=ind; i<tree->size; i++) {
        if(tree->name != NULL)
            tree->name[i] = NULL;
        if(tree->comment != NULL)
            tree->comment[i] = NULL;
    }
    return index;
}

/*remove the subtree pending from n*/
void removeSubtree(int n, TypeTree *tree) {
    int *index = removeSubtreeReturnIndex(n, tree);
    free((void*)index);
}

/*add a new leaf from branch/node n with name n, return the index of the new leaf*/
int addLeaf(int n, char *name, TypeTree *tree) {
    if(tree->size>=tree->sizeBuf)
        reallocTree(tree->sizeBuf+INC_SIZE_TREE, tree);
    tree->node[tree->size].child = n;
    tree->node[tree->size].sibling = tree->node[n].sibling;
    if(tree->time != NULL)
        tree->time[tree->size] = NO_TIME;
    if(tree->parent[n] == NOSUCH) {
        tree->root = tree->size;
    } else {
        int *toSet;
        for(toSet=&(tree->node[tree->parent[n]].child); *toSet!=NOSUCH && *toSet!=n;  toSet = &(tree->node[*toSet].sibling))
            ;
        if(*toSet != n) {
            warning("Issue A while adding new node to  %d (parent %d)\n", n, tree->parent[n]);
            return NOSUCH;
        } else
            *toSet = tree->size;
    }
    tree->parent[tree->size] = tree->parent[n];
    tree->parent[n] = tree->size;
    if(tree->name != NULL)
        tree->name[tree->size] = NULL;
    if(tree->comment != NULL)
        tree->comment[tree->size] = NULL;
    tree->size++;
    if(tree->size>=tree->sizeBuf)
        reallocTree(INC_SIZE_TREE, tree);
    tree->node[tree->size].child = NOSUCH;
    tree->node[tree->size].sibling = NOSUCH;
    tree->node[n].sibling = tree->size;
    if(tree->time != NULL)
        tree->time[tree->size] = NO_TIME;
    tree->parent[tree->size] = tree->size-1;
    if(tree->name != NULL)
        tree->name[tree->size] = name;
    if(tree->comment != NULL)
        tree->comment[tree->size] = NULL;
    tree->size++;
    return tree->size-1;
}



void toBinary(TypeTree *tree) {
    int n;
    for(n=0; n<tree->size; n++) {
        int c = tree->node[n].child;
        if(c != NOSUCH) {
            int d = tree->node[c].sibling;
            if(d != NOSUCH) {
                while(tree->node[d].sibling != NOSUCH) {
                    if(tree->size>=tree->sizeBuf)
                        reallocTree(INC_SIZE_TREE, tree);
                    tree->node[c].sibling = tree->size;
                    tree->node[tree->size].sibling = NOSUCH;
                    tree->node[tree->size].child = d;
                    if(tree->time != NULL)
                        tree->time[tree->size] = NO_TIME;
                    if(tree->parent != NULL) {
                        tree->parent[tree->size] = c;
                        tree->parent[d] = tree->size;
                    }
                    tree->size++;
                    c = d;
                    d = tree->node[d].sibling;
                }
            }
        }
    }
}


/*make node m be child of n*/
void transfer(int m, int n, TypeTree *tree) {
    if(tree == NULL)
        return;
    if(tree->parent == NULL)
        setParent(tree);
    int pm = tree->parent[m], pn;
    if(pm == n)
        return;
    if(pm == NOSUCH)
        return;
    else { //unhang m
        int *toSet;
        for(toSet=&(tree->node[pm].child); *toSet!=NOSUCH && *toSet!=m;  toSet = &(tree->node[*toSet].sibling))
            ;
        if(*toSet != m) {
            warning("Issue A while transfering %d (parent %d) to %d\n", m, pm, n);
            return;
        } else
            *toSet = tree->node[m].sibling;
        if(tree->node[tree->node[pm].child].sibling == NOSUCH) { //pm has to be removed from the tree
            int ppm = tree->parent[pm];
            if(ppm == NOSUCH) { //pm is root
                if(tree->node[pm].child!=NOSUCH) {
                    tree->root = tree->node[pm].child;
                    tree->parent[tree->node[pm].child] = NOSUCH;
                } else {
                    warning("Issue B while transfering %d (parent %d) to %d\n", m, pm, n);
                    return;
                }
            } else {
                int *toSet;
                for(toSet = &(tree->node[ppm].child); *toSet!=NOSUCH && *toSet!=pm; toSet = &(tree->node[*toSet].sibling))
                    ;
                if(*toSet != pm) {
                    warning("Issue C while transfering %d (parent %d, grand %d) to %d\n", m, pm, ppm, n);
                    return;
                } else {
                    tree->node[tree->node[pm].child].sibling = tree->node[pm].sibling;
                    *toSet = tree->node[pm].child;
                    tree->parent[tree->node[pm].child] = ppm;
                }
            }
        } else { //need to add a new node
            if(tree->size>=tree->sizeBuf)
                reallocTree(INC_SIZE_TREE, tree);
            pm = tree->size++;
        }
        tree->node[pm].child = m;
        tree->node[m].sibling = NOSUCH;
        tree->parent[m] = pm;
    }   // from here we have to put pm on the branch n (piece of cake)
    pn = tree->parent[n];
    if(pn == NOSUCH) { //i.e. n = tree->root
        tree->root = pm;
        tree->parent[pm] = NOSUCH;
    } else { // replace n by pm among children of pn
        int *toSet;
        for(toSet=&(tree->node[pn].child); *toSet!=NOSUCH && *toSet!=n; toSet=&(tree->node[*toSet].sibling))
            ;
        if(*toSet != n) {
            warning("Issue D while transfering %d (parent %d) to %d  (parent %d)\n", m, pm, n, pn);
            return;
        } else {
            tree->node[pm].sibling = tree->node[n].sibling;
            *toSet = pm;
            tree->parent[pm] = pn;
        }
    }
    tree->node[pm].child = n;
    tree->node[n].sibling = m;
    tree->parent[n] = pm;
}

/*return 1 if m descends from n*/
int isDescendant(int m, int n, TypeTree *tree) {
    return (getLCA(n, m, tree) == n);
}

/*turn branch lengthes to absolute time rec*/
void bltoabsTimeRec(int n, double offset, TypeTree *tree) {
        int c;
        tree->time[n] = offset;
        for(c=tree->node[n].child; c>=0; c=tree->node[c].sibling)
            bltoabsTimeRec(c, offset+tree->time[c], tree);
}

/*turn branch lengthes to absolute time*/
void bltoabsTime(TypeTree *tree) {
	if(tree->time[tree->root] == NO_TIME)
		tree->time[tree->root] = 0.;
	bltoabsTimeRec(tree->root, tree->time[tree->root], tree);
	if(tree->minTime == NO_TIME)
		tree->minTime = 0.;
}

/*turn branch lengthes to absolute time rec*/
void abstoblTimeRec(int n, double offset, TypeTree *tree) {
        int c;
        for(c=tree->node[n].child; c>=0; c=tree->node[c].sibling)
            abstoblTimeRec(c, tree->time[n], tree);
        tree->time[n] -= offset;
}

/*turn branch lengthes to absolute time*/
void abstoblTime(TypeTree *tree) {
    abstoblTimeRec(tree->root, 0., tree);
}

/*get the smallest time of the tree - makes sense only with absolute times*/
double getMinTime(TypeTree *tree) {
    if(tree->size == 0)
        return 0.;
    return tree->time[tree->root];
}


/*get the greatest time of the tree - makes sense only with absolute times*/
double getMaxTime(TypeTree *tree) {
    int n;
    double max;
    if(tree->size == 0)
        return 0.;
    max = tree->time[0];
    for(n=1; n<tree->size; n++)
        if(tree->time[n]>max)
            max = tree->time[n];
    return max;
}

/*shift all the times/branch from beg*/
void offsetTime(double beg, TypeTree *tree) {
    int i;
    for(i=0; i<tree->size; i++)
        tree->time[i] -= beg;
}

/*return the least common ancestor of n and m*/
int getLCA(int n, int m, TypeTree *tree) {
    int *flag, i, o;
    if(tree->parent == NULL)
        setParent(tree);
    flag = (int*) malloc(tree->size*sizeof(int));
    for(i=0; i<tree->size; i++)
        flag[i] = 0;
    for(o=n; tree->parent[o] >= 0; o = tree->parent[o])
        flag[o] = 1;
    for(o=m; flag[o] == 0 && tree->parent[o] >= 0; o = tree->parent[o])
    ;
    free((void*)flag);
    return o;
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

/*set the table of parents*/
void setParent(TypeTree *tree) {
    int n;
    if(tree == NULL)
        return;
    if(tree->parent == NULL)
        tree->parent = (int*) malloc(tree->sizeBuf*sizeof(int));
    for(n=0; n<tree->sizeBuf; n++)
        tree->parent[n] = NOSUCH;
    for(n=0; n<tree->size; n++) {
        int c;
        for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
            tree->parent[c] = n;
    }
}

/*return the root of tree*/
int getRoot(TypeTree *tree) {
    int n;
    if(tree->size == 0)
        return 0;
    if(tree->parent == NULL)
        setParent(tree);
    for(n=0; n<tree->size && tree->parent[n]>=0; n++);
    return n;
}
/*name (numerote) leaves of tree*/
char **nameBoth(char *prefixIntern, char *prefixLeaf, TypeTree *tree) {
    int n, nLeaves = 1, nInterns = 1, currL = 1, currI = 1, lmax, imax;
    char **name, buffer[200];
    if(tree->size == 0)
        return NULL;
    name = (char**) malloc(tree->size*sizeof(char*));
    for(n=0; n<tree->size; n++)
        if(tree->node[n].child<0)
            nLeaves++;
        else
            nInterns++;

    lmax = (int) floor(log10((double)nLeaves));
    imax = (int) floor(log10((double)nLeaves));
    for(n=0; n<tree->size; n++)
        if(tree->node[n].child<0) {
            char *tmp = buffer;
            int i;
            if(prefixLeaf != NULL)
                tmp += sprintf(tmp, "%s", prefixLeaf);
            for(i=(int) floor(log10((double)currL)); i<lmax; i++)
                tmp += sprintf(tmp, "0");
            tmp += sprintf(tmp, "%d", currL);
            name[n] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
            strcpy(name[n], buffer);
            currL++;
        } else {
            char *tmp = buffer;
            int i;
            if(prefixIntern != NULL)
                tmp += sprintf(tmp, "%s", prefixIntern);
            for(i=(int) floor(log10((double)currI)); i<imax; i++)
                tmp += sprintf(tmp, "0");
            tmp += sprintf(tmp, "%d", currI);
            name[n] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
            strcpy(name[n], buffer);
            currI++;
        }
    return name;
}

/*name (numerote) leaves of tree*/
char **nameLeaves(char *prefix, TypeTree *tree) {
    int n, nLeaves = 1, curr = 1, lmax;
    char **name, buffer[200];
    if(tree->size == 0)
        return NULL;
    name = (char**) malloc(tree->sizeBuf*sizeof(char*));
    for(n=0; n<tree->size; n++)
        if(tree->node[n].child<0)
            nLeaves++;
    lmax = (int) floor(log10((double)nLeaves));
    for(n=0; n<tree->size; n++)
        if(tree->node[n].child<0) {
            char *tmp = buffer;
            int i;
            if(prefix != NULL)
                tmp += sprintf(tmp, "%s", prefix);
            for(i=(int) floor(log10((double)curr)); i<lmax; i++)
                tmp += sprintf(tmp, "0");
            tmp += sprintf(tmp, "%d", curr);
            name[n] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
            strcpy(name[n], buffer);
            curr++;
        } else
            name[n] = NULL;
    for(n=tree->size; n<tree->sizeBuf; n++)
        name[n] = NULL;
    return name;
}

/*return a standard feature
TypeStandardFeature *getBasicStandardFeature(TypeTree *tree) {
    TypeStandardFeature *res;
    res = (TypeStandardFeature*) malloc(sizeof(TypeStandardFeature));
    res->name = nameBoth("Node_", "Leaf_", tree);
    res->comment = NULL;
    res->info = NULL;
    return res;
}*/

/*get the status from comments*/
char *getSpecy(char *str) {
    int i, ind = 0;
    char tmp[BASIC_TMP_SIZE], *specy;

    if(str == NULL || strlen(str) < 2)
        return NULL;
    for(i=2; str[i] != '\0' && (str[i-2]!=':' || str[i-1]!='S' || str[i]!='='); i++);
    i++;
    for(; str[i] != '\0' && issep(str[i]); i++);
    for(; i<BASIC_TMP_SIZE-1 && str[i] != '\0' && str[i] != ':' && !issep(str[i]); i++)
        tmp[ind++] = str[i];
    tmp[ind++] = '\0';
    if(ind>1) {
        specy = (char*) malloc((strlen(tmp)+1)*sizeof(char));
        strcpy(specy, tmp);
    } else
        specy = NULL;
    return specy;
}

char **getNamestoSpecies(char **comment, TypeTree *tree) {
    int n;
    char **name;
    name = (char**) malloc(tree->size*sizeof(char*));
    for(n=0; n<tree->size; n++) {
        name[n] = getSpecy(comment[n]);
    }
    return name;
}

char **getNamestoSpeciesLeaves(char **comment, TypeTree *tree) {
    int n;
    char **name;
    name = (char**) malloc(tree->size*sizeof(char*));
    for(n=0; n<tree->size; n++) {
        if(tree->node[n].child < 0)
            name[n] = getSpecy(comment[n]);
        else
            name[n] = NULL;
    }
    return name;
}

/*returns the number of contemporary lineages (i.e. living at maxTime) of "tree"*/
int countContemp(TypeTree *tree) {
    int n, count = 0;
    if(tree == NULL)
           return 0;
    for(n=0; n<tree->size; n++)
        if(tree->time[n]>=tree->maxTime)
            count++;
    return count;
}

/*returns the number of leaves of "tree"*/
int countLeaves(TypeTree *tree) {
    int n, count = 0;
    if(tree == NULL)
           return 0;
    for(n=0; n<tree->size; n++)
        if(tree->node[n].child<0)
            count++;
    return count;
}

/*returns all the leaves of "tree"*/
int *getAllLeaves(TypeTree *tree) {
  int n, count = 0;
  if(tree == NULL) return NULL;
  int *allLeaves;
  allLeaves = (int*) malloc(countLeaves(tree));
  for(n=0; n<tree->size; n++) {
    if(tree->node[n].child<0) {
      allLeaves[count] = n;
      count++;
    }
  }
  return allLeaves;
}

/*return the numbre of childre of node n*/
int getNumberChildren(int n, TypeTree *tree) {
    int c, count = 0;
    for(c=tree->node[n].child; c != -1; c=tree->node[c].sibling)
        count ++;
    return count;
}

/*returns the number of leaves of the n subtree of "tree"*/
int countSubLeaves(int n, TypeTree *tree) {
    int c, count = 0;
    if(tree->node[n].child == -1)
        return 1;
    for(c=tree->node[n].child; c != -1; c=tree->node[c].sibling)
        count += countSubLeaves(c, tree);
    return count;
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
    res->maxTime = tree->maxTime;
    res->minTime = tree->minTime;
    res->maxTimeInt = tree->maxTimeInt;
    res->minTimeInt = tree->minTimeInt;
    return res;
}

void reallocTree(int size, TypeTree *tree) {
    tree->node = (TypeNode*) realloc((void*) tree->node, size*sizeof(TypeNode));
    if(tree->time != NULL) {
        int i;
        tree->time = (double*) realloc((void*) tree->time, size*sizeof(double));
        for(i=tree->sizeBuf; i<size; i++)
            tree->time[i] = NO_TIME;
    }
    if(tree->parent != NULL) {
        int i;
        tree->parent = (int*) realloc((void*) tree->parent, size*sizeof(int));
        for(i=tree->sizeBuf; i<size; i++)
            tree->parent[i] = NOSUCH;
    }
    if(tree->name != NULL) {
        int i;
        tree->name = (char**) realloc((void*) tree->name, size*sizeof(char*));
        for(i=tree->sizeBuf; i<size; i++)
            tree->name[i] = NULL;
    }
    if(tree->comment != NULL) {
        int i;
        tree->comment = (char**) realloc((void*) tree->comment, size*sizeof(char*));
        for(i=tree->sizeBuf; i<size; i++)
            tree->comment[i] = NULL;
    }
    tree->sizeBuf = size;
}

/*allocate a string tab and set its entries to NULL*/
char **newStringTab(int size) {
    int i;
    char **res;
    res = (char**) malloc(size*sizeof(char*));
    for(i=0; i<size; i++)
        res[i] = NULL;
    return res;
}

/*allocate and initialize a new_feat tree*/
TypeTree *newTree(int sizeBuf) {
    TypeTree *tree;
    tree = (TypeTree*) malloc(sizeof(TypeTree));
    tree->sizeBuf = sizeBuf;
    if(tree->sizeBuf > 0) {
        int i;
        tree->node = (TypeNode*) malloc(tree->sizeBuf*sizeof(TypeNode));
        tree->time = (double*) malloc(sizeBuf*sizeof(double));
        for(i=0; i<tree->sizeBuf; i++)
            tree->time[i] = NO_TIME;
    } else {
        tree->node = NULL;
        tree->time = NULL;
    }
    tree->parent = NULL;
    tree->info = NULL;
    tree->name = NULL;
    tree->comment = NULL;
    tree->size = 0;
    tree->minTime = NO_TIME;
    tree->minTimeInt.inf = NO_TIME;
    tree->minTimeInt.sup = NO_TIME;
    tree->maxTime = NO_TIME;
    tree->maxTimeInt.inf = NO_TIME;
    tree->maxTimeInt.sup = NO_TIME;
    return tree;
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

/*set entry n of dest to entry m of src*/
void setTree(int n, TypeTree *dest, int m, TypeTree *src) {
    dest->node[n] = src->node[m];
    if(src->time != NULL && dest->time != NULL)
        dest->time[n] = src->time[m];
}



/*initialize fields of node to standard values*/
void initNode(TypeNode *node) {
    node->child = -1;
    node->sibling = -1;
}



/*read tree in newick format*/
TypeTree *readTreeFromName(char *filename) {
    FILE *f;
    if((f = fopen(filename, "r")))
        return readTree(f);
    else
        return NULL;
}

/*read tree in newick format*/
TypeTree *readTree(FILE *f) {
    char c;
    TypeTree *tree;
    tree = newTree(INC_SIZE);
    tree->name = newStringTab(INC_SIZE);
    tree->comment = newStringTab(INC_SIZE);
    skipSeparator(f);
    c = fgetc(f);
    if(c == '(' ) {
        tree->root = readNode(f, tree);
    }
    skipSeparator(f);
    c = fgetc(f);
    if(c != ';') {
        error("error not ending by ; - current caracter '%c'\n", c);
    }
    for(c=0; c<tree->size && tree->name == NULL; c++)
        ;
    if(c == tree->size) {
        free((void*)tree->name);
        tree->name = NULL;
    }
    for(c=0; c<tree->size && tree->comment == NULL; c++)
        ;
    if(c == tree->size) {
        free((void*)tree->comment);
        tree->comment = NULL;
    }
    return tree;
}

#define INC_SIZE_TAB 10
/*read tree in newick format*/
TypeTree **readTrees(FILE *f) {
    char c;
    int size, sizeBuf;
    TypeTree **tree;
     size = 0;
    sizeBuf = INC_SIZE_TAB;
    tree = malloc(sizeBuf*sizeof(TypeTree*));
    skipSeparator(f);
    do {
        if(size>=sizeBuf) {
            sizeBuf += INC_SIZE_TAB;
            tree = realloc((void*) tree, sizeBuf*sizeof(TypeTree*));
         }
        tree[size] = newTree(INC_SIZE);
        tree[size]->name = newStringTab(INC_SIZE);
        tree[size]->comment = newStringTab(INC_SIZE);
         c = fgetc(f);
        if(c == '(' ) {
            tree[size]->root = readNode(f, tree[size]);
        }
        skipSeparator(f);
        c = fgetc(f);
        if(c != ';') {
            freeTree(tree[size]);
         } else {
            size++;
        }
        skipSeparator(f);
    } while(c != EOF);
    tree = realloc((void*) tree, (size+1)*sizeof(TypeTree*));
    tree[size] = NULL;
    return tree;
}
/*read a node in newick format*/
int readNode(FILE *f, TypeTree *tree) {
    int current, child;
    char c;

    if(tree->size >= tree->sizeBuf)
        reallocTree(tree->sizeBuf+INC_SIZE, tree);
    current = tree->size++;
    initNode(&(tree->node[current]));
    do {
        skipSeparator(f);
        c = fgetc(f);
        if(c == '(') {
            child = readNode(f, tree);
            tree->node[child].sibling = tree->node[current].child;
            tree->node[current].child = child;
        } else {
            ungetc(c, f);
            child = readLeaf(f, tree);
            tree->node[child].sibling = tree->node[current].child;
            tree->node[current].child = child;
        }
        skipSeparator(f);
        c = fgetc(f);
    } while(c == ',');
    if(c == ')') {
        tree->name[current] = readName(f);
        tree->time[current] = readTime(f);
        tree->comment[current] = readComment(f);
    } else {
        error("unclosed parenthesis - current caracter '%c'\n", c);
    }
    return current;
}

/*read a leaf in newick format*/
int readLeaf(FILE *f, TypeTree *tree) {
    int current;
    if(tree->size >= tree->sizeBuf)
        reallocTree(tree->sizeBuf+INC_SIZE, tree);
    current = tree->size++;
    initNode(&(tree->node[current]));
    tree->node[current].child = -1;
    tree->name[current] = readName(f);
    tree->time[current] = readTime(f);
    tree->comment[current] = readComment(f);
    return current;
}

/*read the name*/
char *readName(FILE *f) {
    char c, *tmp;
    int i;
    tmp = (char*) malloc((BASIC_TMP_SIZE+1)*sizeof(char));
    skipSeparator(f);
    c = fgetc(f);
    if(c == '\'' || c == '"') {
        c = fgetc(f);
        for(i=0; i<BASIC_TMP_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
            tmp[i] = c;
            c = fgetc(f);
        }
        if(c == '\'' || c == '"')
            c = fgetc(f);
        else
            error("Missing closing \' or \"...");
    } else {
        for(i=0; i<BASIC_TMP_SIZE && c !=EOF && c != ':' && c != '(' && c != ')'  && c != '[' && c != ',' && c != ';'; i++) {
            tmp[i] = c;
            c = fgetc(f);
        }
    }
    if(i == BASIC_TMP_SIZE)
        error("Name too much long...");
    ungetc(c, f);
    tmp[i++] = '\0';
    removeSpaces(tmp);
    if(strlen(tmp)>0)
        return (char *) realloc((void*) tmp, i*sizeof(char));
    else {
        free((void*) tmp);
        return NULL;
    }
}

/*read the comment*/
char *readComment(FILE *f) {
    char c, *tmp;
    int i=0;
    tmp = (char*) malloc((BASIC_TMP_SIZE+1)*sizeof(char));
    skipSeparator(f);
    c = fgetc(f);
    if(c=='[') {
        c = fgetc(f);
        for(i=0; i<BASIC_TMP_SIZE && c !=EOF && c != ']'; i++) {
            tmp[i] = c;
            c = fgetc(f);
        }
        if(i == BASIC_TMP_SIZE)
            error("comment too much long...");
    } else
        ungetc(c, f);
    tmp[i++] = '\0';
    if(i>1)
        return (char *) realloc((void*) tmp, i*sizeof(char));
    else {
        free((void*) tmp);
        return NULL;
    }
}

/*read the time or branch length*/
double readTime(FILE *f) {
    char tmp[BASIC_TMP_SIZE+1], c;
    int i;
    double time = NO_TIME;
    skipSeparator(f);
    c = fgetc(f);
    if(c==':') {
        c = fgetc(f);
        skipSeparator(f);
        for(i=0; i<BASIC_TMP_SIZE && c !=EOF && c != '(' && c != ')'&& c != ',' && c != '[' && c != ';'; i++) {
            tmp[i] = c;
            c = fgetc(f);
        }
        ungetc(c,f);
        tmp[i] = '\0';
        time = atof(tmp);
    } else
        ungetc(c, f);
    return time;
}


/*compare two reorder type*/
int compareTmpReorder(const void *a, const void *b) {
    if(((TypeTmpReorder*)a)->name == NULL) {
        if(((TypeTmpReorder*)b)->name == NULL)
            return 0;
        else
            return -1;
    } else {
        if(((TypeTmpReorder*)b)->name == NULL)
            return 1;
        else
            return strcmp(((TypeTmpReorder*)a)->name,((TypeTmpReorder*)b)->name);
    }
}

void reorderTreeRec(char **tmp, char **name, int n, TypeTree *tree) {
    int c, nchild = 0;
    for(c=tree->node[n].child; c>=0; c=tree->node[c].sibling) {
        nchild++;
        reorderTreeRec(tmp, name, c, tree);
    }
    if(nchild>0) {
        int i;
        TypeTmpReorder *tab;
        tab = (TypeTmpReorder*) malloc(nchild*sizeof(TypeTmpReorder));
        nchild = 0;
        for(c=tree->node[n].child; c>=0; c=tree->node[c].sibling) {
            tab[nchild].name = tmp[c];
            tab[nchild].index = c;
            nchild++;
        }
        qsort(tab, nchild, sizeof(TypeTmpReorder), compareTmpReorder);
        tree->node[n].child = tab[0].index;
        for(i=1; i<nchild; i++)
            tree->node[tab[i-1].index].sibling = tab[i].index;
        tree->node[tab[nchild-1].index].sibling = -1;
        free((void*)tab);
    }
    if(nchild>0)
        tmp[n] = tmp[tree->node[n].child];
    else {
        if(name && name[n])
            tmp[n] = name[n];
        else
            tmp[n] = NULL;
    }
/*
 * an alternative order taking into account name od internal nodes
 * 	if(tree->name[n] != NULL)
        name[n] = tree->name[n];
    else {
        if(nchild>0)
            name[n] = name[tree->node[n].child];
        else
            name[n] = NULL;
    }
*/}

/*reorder all children in a canonical (provided at least leaves are named)*/
void reorderTree(char **name, TypeTree *tree) {
    char **tmp;
    if(tree == NULL || name == NULL)
        return;
    tmp = (char**) malloc(tree->size*sizeof(char*));
    reorderTreeRec(tmp, name, tree->root, tree);
    free((void*)tmp);
}

/*reorder all children in a canonical (provided at least leaves are named)*/
TypeTree *reorderTreeCpy(char **name, TypeTree *tree) {
    TypeTree *reor;
    char **tmp;
    if(tree == NULL || name == NULL)
        return NULL;
    reor = cpyTree(tree);
    tmp = (char**) malloc(reor->size*sizeof(char*));
    reorderTreeRec(tmp, name, reor->root, reor);
    free((void*)tmp);
    return reor;
}

/*return a tree where leaves appear in lexicographic order*/
TypeLexiTree *getDictTree(char **name, int size) {
    int i;
    TypeLexiTree *dict;
    dict = newLexiTree();
    for(i=0; i<size; i++)
        if(name && name[i]) {
            if(addWordLexi(name[i], i, dict)>=0)
                warning("Warning! duplicate identifier '%s'\n", name[i]);
        }
    return dict;
}


char *nameInternalNodesRec(int n, char **name, TypeTree *tree) {
    int c;
    char *tmp, *min;
    if(tree->node[n].child < 0)
        return name[n];
    tmp = (char*) malloc(BASIC_TMP_SIZE*sizeof(char));
    min = nameInternalNodesRec(tree->node[n].child, name, tree);
    if(min)
        strcpy(tmp, min);
    for(c=tree->node[tree->node[n].child].sibling; c>=0; c=tree->node[c].sibling) {
        char *childName;
        childName = nameInternalNodesRec(c, name, tree);
        strcat(tmp, "_*_");
        if(childName) {
            strcat(tmp, childName);
            if(min == NULL || strcmp(childName, min) < 0)
                min = childName;
        }
    }
    tmp = (char*) realloc((void*)tmp, (strlen(tmp)+1)*sizeof(char));
    if(name[n] != NULL)
        free((void*) name);
    name[n] = tmp;
    return min;
}

void nameInternalNodes(char **name, TypeTree *tree) {
    nameInternalNodesRec(tree->root, name, tree);
}

TypeTimeTab *readTimeTab(FILE *f) {
    int sizeBuf;
    char c;
    TypeTimeTab *res;

    res = (TypeTimeTab *) malloc(sizeof(TypeTimeTab));
    sizeBuf = BASIC_INC_BUFFER;
    res->name1 = (char**) malloc(sizeBuf*sizeof(char*));
    res->name2 = (char**) malloc(sizeBuf*sizeof(char*));
    res->time = (double*) malloc(sizeBuf*sizeof(double));
    res->size = 0;
    do {
        for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
        if(c != EOF) {
            char *tmp;
            int i;
            if(res->size >= sizeBuf) {
                sizeBuf += BASIC_INC_BUFFER;
                res->name1 = (char**) realloc((void*) res->name1, sizeBuf*sizeof(char*));
                res->name2 = (char**) realloc((void*) res->name2, sizeBuf*sizeof(char*));
                res->time = (double*) realloc((void*) res->name2, sizeBuf*sizeof(double));
            }
            tmp = (char*) malloc((BASIC_TMP_SIZE+1)*sizeof(char));
            if(c == '\'' || c == '"') {
                c = fgetc(f);
                for(i=0; i<BASIC_TMP_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
                    tmp[i] = c;
                    c = fgetc(f);
                }
                if(c == '\'' || c == '"')
                    c = fgetc(f);
            } else {
                for(i=0; i<BASIC_TMP_SIZE && c !=EOF && !issep(c); i++) {
                    tmp[i] = c;
                    c = fgetc(f);
                }
            }
            if(i == BASIC_TMP_SIZE)
                error("Name too much long...");
            tmp[i++] = '\0';
            res->name1[res->size] = (char *) realloc((void*) tmp, i*sizeof(char));
            for(; c != EOF && issep(c); c = fgetc(f));
            tmp = (char*) malloc((BASIC_TMP_SIZE+1)*sizeof(char));
            if(c == '\'' || c == '"') {
                c = fgetc(f);
                for(i=0; i<BASIC_TMP_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
                    tmp[i] = c;
                    c = fgetc(f);
                }
                if(c == '\'' || c == '"')
                    c = fgetc(f);
            } else {
                for(i=0; i<BASIC_TMP_SIZE && c !=EOF && !issep(c); i++) {
                    tmp[i] = c;
                    c = fgetc(f);
                }
            }
            if(i == BASIC_TMP_SIZE)
                error("Name too much long...");
            tmp[i++] = '\0';
            res->name2[res->size] = (char *) realloc((void*) tmp, i*sizeof(char));
            for(; c != EOF && issep(c); c = fgetc(f));
            tmp = (char*) malloc((BASIC_TMP_SIZE+1)*sizeof(char));
            for(i=0; c != EOF && !issepline(c) && i<BASIC_TMP_SIZE; i++) {
                tmp[i] = c;
                c = fgetc(f);
            }
            if(i == BASIC_TMP_SIZE)
                error("Time too much long...");
            tmp[i++] = '\0';
            res->time[res->size] = atof(tmp);
            free((void*)tmp);
            res->size++;
        }
    } while(c != EOF);
    return res;
}

/*read a table <leaf1>\t<leaf2>\t<time> and set the time of the lca of <leaf1> and <leaf2> to <time>*/
void setTime(FILE *fi, char **name, TypeTree *tree) {
	TypeLexiTree *dict;
	TypeTimeTab *list;
	int i, n, m;
	dict = getDictTree(name, tree->size);
	list = readTimeTab(fi);
	if(tree->parent == NULL)
		setParent(tree);
	for(i=0; i<list->size; i++)
		if(((n = findWordLexi(list->name1[i], dict))>=0) && ((m = findWordLexi(list->name2[i], dict))>=0))
			tree->time[getLCA(n, m, tree)] = list->time[i];
		else
			warning("Problem with reading line %s\t%s\n", list->name1[i], list->name2[i]);
	free((void*)list->name1);
	free((void*)list->name2);
	free((void*)list->time);
	free((void*)list);
	freeLexiTree(dict);
}



/*prune "tree" to keep only leaves in dict*/
TypeTree *pruneLeavesFromDict(TypeTree *tree, TypeLexiTree *dict) {
	TypeTree *resT;
	int n, *new, *parent;
	if(tree->name != NULL) {
		resT = newTree(tree->size);
		resT->name = (char**) malloc(tree->size*sizeof(char*));;
		resT->maxTime = tree->maxTime;
		resT->minTime = tree->minTime;
	} else {
		resT = newTree(0);
		return resT;
	}
	resT->size = 0;
	new = (int*) malloc(tree->size*sizeof(int));
	for(n=0; n<tree->size; n++)
		new[n] = NOSUCH;
	parent = getParent(tree);
	for(n=0; n<tree->size; n++) {
		if(tree->node[n].child == NOSUCH && tree->name[n] != NULL && findWordLexi(tree->name[n], dict) != -1) {
			int m, prec;
			new[n] = resT->size++;
			prec = new[n];
			resT->node[new[n]].child = NOSUCH;
			resT->node[new[n]].sibling = NOSUCH;
			resT->time[new[n]] = tree->time[n];
			resT->name[new[n]] = strdpl(tree->name[n]);
			for(m=parent[n]; m != NOSUCH && new[m] == NOSUCH; m=parent[m]) {
				new[m] = resT->size++;
				resT->node[new[m]].child = prec;
				prec = new[m];
				resT->node[new[m]].sibling = NOSUCH;
				if(tree->time != NULL)
					resT->time[new[m]] = tree->time[m];
				resT->name[new[m]] = strdpl(tree->name[m]);
			}
			if(m != NOSUCH) {
				resT->node[prec].sibling = resT->node[new[m]].child;
				resT->node[new[m]].child = prec;
			}
		}
	}
	resT->sizeBuf = resT->size;
	resT->name = (char**) realloc(resT->name, resT->size*sizeof(char*));
	resT->node = (TypeNode*) realloc(resT->node, resT->size*sizeof(TypeNode));
	if(resT->time)
		resT->time = (double*) realloc(resT->time, resT->size*sizeof(double));
	free((void*)new);
	free((void*)parent);
	parent = getParent(resT);
	for(n=0; n<resT->size && parent[n]>=0; n++);
	free((void*)parent);
	resT->root = n;
	return resT;
}

/*prune "tree" to that can be observed from contemporary lineages only*/
TypeTree *pruneContemp(TypeTree *tree) {
    TypeTree *resT;
    int n, *new_feat;
    if(tree->time == NULL)
        return NULL;
    resT = newTree(tree->sizeBuf);
    resT->maxTime = tree->maxTime;
    resT->minTime = tree->minTime;
    resT->size = 0;
    if(tree->size == 0)
        return resT;
    if(tree->name)
        resT->name = newStringTab(resT->sizeBuf);
    if(tree->comment)
        resT->comment = newStringTab(resT->sizeBuf);
    new_feat = (int*) malloc(tree->size*sizeof(int));
    for(n=0; n<tree->size; n++)
        new_feat[n] = NOSUCH;
    if(tree->parent == NULL)
        setParent(tree);
    for(n=0; n<tree->size; n++) {
        if(tree->time[n]>=tree->maxTime) {
            int m;
            new_feat[n] = resT->size;
            resT->node[resT->size].child = NOSUCH;
            resT->node[resT->size].sibling = NOSUCH;
            resT->time[resT->size] = tree->time[n];
			if(tree->name) {
				if(tree->name[n] != NULL)
					resT->name[resT->size] = strdpl(tree->name[n]);
				else
					resT->name[resT->size] = NULL;
			}
			if(tree->comment) {
				if(tree->comment[n] != NULL) 
					resT->comment[resT->size] = strdpl(tree->comment[n]);
				else
					resT->comment[resT->size] = NULL;
			}
           resT->time[resT->size] = tree->time[n];
            resT->size++;
            for(m=tree->parent[n]; m>=0 && new_feat[m]==NOSUCH; m=tree->parent[m]) {
                new_feat[m] = resT->size;
                resT->node[resT->size].child = resT->size-1;
                resT->node[resT->size].sibling = NOSUCH;
                resT->time[resT->size] = tree->time[m];
				if(tree->name) {
					if(tree->name[n] != NULL)
						resT->name[resT->size] = strdpl(tree->name[n]);
					else
						resT->name[resT->size] = NULL;
				}
				if(tree->comment) {
					if(tree->comment[n] != NULL) 
						resT->comment[resT->size] = strdpl(tree->comment[n]);
					else
						resT->comment[resT->size] = NULL;
				}
                resT->size++;
            }
            if(m>=0) {
                resT->time[new_feat[m]] = tree->time[m];
                if(resT->node[new_feat[m]].child>=0) {
                    resT->node[resT->node[new_feat[m]].child].sibling = resT->size-1;
                } else {
                    resT->node[new_feat[m]].child = resT->size-1;
                }
            }
        }
    }
    free((void*)new_feat);
    if(resT->parent == NULL)
        setParent(resT);
    for(n=0; n<resT->size && resT->parent[n]!=NOSUCH; n++);
    resT->root = n;
    return resT;
}



int iterateBinary(int n, TypeTree *resT, TypeTree *tree) {
    int m;
    for(m=n; tree->node[m].child>=0 && tree->node[tree->node[m].child].sibling<0; m=tree->node[m].child);
    if(tree->node[m].child>=0) {
        int c1, c2, c;
        c1 = iterateBinary(tree->node[m].child, resT, tree);
        c2 = c1;
        for(c=tree->node[tree->node[m].child].sibling; c >= 0; c = tree->node[c].sibling) {
            resT->node[c2].sibling = iterateBinary(c, resT, tree);
            c2 = resT->node[c2].sibling;
        }
        resT->node[c2].sibling = -1;
        resT->node[resT->size].child = c1;
    } else {
        resT->node[resT->size].child = -1;
    }
    if(tree->time)
        resT->time[resT->size] = tree->time[m];
    if(tree->name)
        resT->name[resT->size] = strdpl(tree->name[m]);
    if(tree->comment)
        resT->comment[resT->size] = strdpl(tree->comment[m]);
    resT->size++;
    return resT->size-1;
}

TypeTree *fixBinary(TypeTree *tree) {
    TypeTree *resT;
    int n;
    resT = cpyTree(tree);
    if(resT->name != NULL)
        for(n=0; n<resT->sizeBuf; n++) {
            if(resT->name[n] != NULL)
                free((void*)resT->name[n]);
            resT->name[n] = NULL;
        }
    if(resT->comment != NULL)
        for(n=0; n<resT->sizeBuf; n++) {
            if(resT->comment[n] != NULL)
                free((void*)resT->comment[n]);
            resT->comment[n] = NULL;
        }
    resT->size = 0;
    iterateBinary(tree->root, resT, tree);
    if(resT->parent == NULL)
		resT->parent = (int*) malloc(resT->sizeBuf*sizeof(int));
	setParent(resT);
    for(n=0; n<resT->size && resT->parent[n]!=NOSUCH; n++);
    resT->root = n;
    return resT;
}



/*replace space by underscore*/
void fixNameUnderscore(char **name, int size) {
    int n;
    for(n=0; n<size; n++)
        if(name && name[n]) {
            replaceChar(name[n],'_', ' ');
            replaceChar(name[n],'\'', ' ');
            replaceChar(name[n],'"', ' ');
        }
}


/*print tree in newick format*/
void fprintSubtreeNewick(FILE *f, int n, TypeTree *tree) {
    if(tree->size<=0 || n>=tree->size || n<0)
        return;
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
    fprintf(f, ";\n");
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

/*print the suffix tree*/
void printTreeDebug(FILE *f, int n, TypeTree *tree, char **name) {
    printNodeDebug(f, n, 0, tree, name);
}


/*print recursively the suffix node*/
void printNodeDebug(FILE *f, int s, int depth, TypeTree *tree, char **name) {
    int tmp;

    if(s < 0)
        return;
    if(depth>=0) {
        int d = depth;
        /*Print the branches coming from higher nodes.*/
        for(d=0; d<depth; d++)
            fprintf(f, "|");
        fprintf(f, "+");
        if(name != NULL && name[s] != NULL)
            fprintf(f, "%s", name[s]);
        else
            fprintf(f, "%d", s);
        if(tree->time[s] != NEG_INFTY)
            fprintf(f, " %.2lf\n", tree->time[s]);
        else
            fprintf(f, " -\n");
    }
    for(tmp=tree->node[s].child; tmp>=0; tmp=tree->node[tmp].sibling)
        printNodeDebug(f, tmp, depth+1, tree, name);
}


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
void fprintTree(FILE *f, TypeTree *tree, TypeDisplayName display) {
	if(tree->size<=0)
		return;
	if(tree->node[tree->root].child >= 0) {
		int tmp = tree->node[tree->root].child;
		fprintf(f, "(");
		fprintNode(f, tmp, tree, display);
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
			fprintf(f, ", ");
			fprintNode(f, tmp, tree, display);
		}
		fprintf(f, ")");
	}
	fprintIdentTimeComment(f, tree->root, tree, display);
	fprintf(f, ";\n");
}


/*print node in newick format*/	
void fprintNode(FILE *f, int n, TypeTree *tree, TypeDisplayName display) {
	if(tree->node[n].child >= 0) {
		int tmp = tree->node[n].child;
		fprintf(f, "(");
		fprintNode(f, tmp, tree, display);
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
			fprintf(f, ", ");
			fprintNode(f, tmp, tree, display);
		}
		fprintf(f, ")");
	}
	fprintIdentTimeComment(f, n, tree, display);
}


/*print tree in PsTricks format*/	
void fprintTreePst(FILE *f, TypeTree *tree, TypeDisplayName display) {
	if(tree->size<=0)
		return;
	if(tree->node[tree->root].child >= 0) {
		int tmp = tree->node[tree->root].child;
		fprintf(f, "\\begin{pspicture}[showgrid=false](-8,-8)(8,8)$\\pstree[treemode=R]{\\TR{");
		fprintIdentTimeComment(f, tree->root, tree, display);
		fprintf(f, "}}{");
		fprintNodePst(f, tmp, tree, display);
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
			fprintf(f, "\n");
			fprintNodePst(f, tmp, tree, display);
		}
		fprintf(f, "}$\\end{pspicture}");
	}
}


/*print node in newick format*/	
void fprintNodePst(FILE *f, int n, TypeTree *tree, TypeDisplayName display) {
	if(tree->node[n].child >= 0) {
		int tmp = tree->node[n].child;
		fprintf(f, "\\pstree[treemode=R]{\\TR{");
		fprintIdentTimeComment(f, n, tree, display);
		fprintf(f, "}}{");
		fprintNodePst(f, tmp, tree, display);
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
			fprintf(f, "\n ");
			fprintNodePst(f, tmp, tree, display);
		}
		fprintf(f, "}");
	} else {
		fprintf(f, "\\TR{");
		fprintIdentTimeComment(f, n, tree, display);
		fprintf(f, "}");
	}
}


/*print tree in debug mode*/
void fprintTreeX(FILE *f, TypeTree *tree) {
    int n;
    fprintf(f, "root %d\n", tree->root);
    if(tree->size==0)
        warning("Empty tree\n");
    else
        for(n=0; n<tree->size; n++) {
            fprintf(f, "%d", n);
             if(tree->name != NULL) {
                    fprintf(f, " (%s) ", tree->name[n]);
            }
           if(tree->time != NULL) {
                if(tree->time[n] == NO_TIME)
                    fprintf(f, ":? ");
                else
                    fprintf(f, ":%.2lf ", tree->time[n]);
            }
            if(tree->node[n].child >= 0) {
                int c;
                fprintf(f, "-> %d", tree->node[n].child);
                for(c=tree->node[tree->node[n].child].sibling; c!=NOSUCH; c=tree->node[c].sibling)
                    fprintf(f, " - %d", c);
            }
            fprintf(f, "\n");
         }
}


void reindexTree(TypeTree *tree, int *index) {
	int n;
	TypeNode *oldNode;
	double *oldTime;
	char **oldName, **oldComment;

	oldNode = tree->node;
	tree->node = (TypeNode*) malloc(tree->sizeBuf*sizeof(TypeNode));
	oldTime = tree->time;
	if(tree->time != NULL)
		tree->time = (double*) malloc(tree->sizeBuf*sizeof(double));
	oldName = tree->name;
	if(tree->name != NULL) {
		tree->name = (char**) malloc(tree->sizeBuf*sizeof(char*));
		for(n=tree->size; n<tree->sizeBuf ; n++)
			tree->name[n] = NULL;
	}
	oldComment = tree->comment;
	if(tree->comment != NULL) {
		tree->comment = (char**) malloc(tree->sizeBuf*sizeof(char*));
		for(n=tree->size; n<tree->sizeBuf ; n++)
			tree->comment[n] = NULL;
	}
	tree->root = index[tree->root];    
	for(n=0; n<tree->size ; n++) {
		tree->node[index[n]].child = (oldNode[n].child != NOSUCH) ? index[oldNode[n].child] : NOSUCH;
		tree->node[index[n]].sibling = (oldNode[n].sibling != NOSUCH) ? index[oldNode[n].sibling] : NOSUCH;
		if(tree->time != NULL)
			tree->time[index[n]] = oldTime[n];
		if(tree->name != NULL)
			tree->name[index[n]] = oldName[n];
		if(tree->comment != NULL)
			tree->comment[index[n]] = oldComment[n];
	}
	free((void*)oldNode);
	if(oldTime != NULL)
		free((void*)oldTime);
	if(oldName != NULL)
		free((void*)oldName);
	if(oldComment != NULL)
		free((void*)oldComment);
}

void reindexLeavesFirst(TypeTree *tree) {
	int *index, le, in, n;
	index = (int*) malloc(tree->size*sizeof(int));
	le = 0; 
	in = tree->size-1;
	for(n=0; n<tree->size ; n++)
		if(tree->node[n].child == NOSUCH)
			index[n] = le++;  
		else
			index[n] = in--;
	reindexTree(tree, index);
	free((void*)index);
}

void fillStandardOrderTip(int n, TypeTree *tree, int *index, int *cur) {
	int c;
	if(tree->node[n].child != NOSUCH) {
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling)
			fillStandardOrderTip(c, tree, index, cur);
	 } else {

		index[n] = (*cur)++;
	}
}

void fillStandardOrderInternal(int n, TypeTree *tree, int *index, int *cur) {
	int c;
	if(tree->node[n].child != NOSUCH) {
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling)
			fillStandardOrderInternal(c, tree, index, cur);
		index[n] = (*cur)++;
	}
}



/*the order returned is such that the tips have the lowest orders, the tips of a subtree have contiguous orders, if a node is ancestor of antother then its order is greater*/
int *getStandardOrder(TypeTree *tree) {
	int cur=0, *index;
	index = (int*) malloc(tree->size*sizeof(int));
	fillStandardOrderTip(tree->root, tree, index, &cur);
	fillStandardOrderInternal(tree->root, tree, index, &cur);
	return index;
}

double *reindexCpyDoubleTable(int *index, double *table, int size) {
	int i;
	double *res;
	res = (double*) malloc(size*sizeof(double));
	for(i=0; i<size; i++)
		res[index[i]] = table[i];
	return res;
}

TypeTree *reindexCpyTree(int *index, TypeTree *tree) {
	TypeTree *res;
	int n;
	res = newTree(tree->size);
	if(tree->name != NULL)
		res->name = (char**) malloc(tree->size*sizeof(char*));
	if(tree->comment != NULL)
		res->comment = (char**) malloc(tree->size*sizeof(char*));
	for(n=0; n<tree->size; n++) {
		if(tree->node[n].child != NOSUCH)
			res->node[index[n]].child = index[tree->node[n].child];
		else
			res->node[index[n]].child = NOSUCH;
		if(tree->node[n].sibling != NOSUCH)
			res->node[index[n]].sibling = index[tree->node[n].sibling];
		else
			res->node[index[n]].sibling = NOSUCH;
		res->time[index[n]] = tree->time[n];
		if(tree->name != NULL) {
			if(tree->name[n] != NULL) {
				res->name[index[n]] = (char*) malloc((strlen(tree->name[n])+1)*sizeof(char));
				strcpy(res->name[index[n]], tree->name[n]);
			} else
				res->name[index[n]] = NULL;
		}
		if(tree->comment != NULL) {
			if(tree->comment[n] != NULL) {
				res->comment[index[n]] = (char*) malloc((strlen(tree->comment[n])+1)*sizeof(char));
				strcpy(res->comment[index[n]], tree->comment[n]);
			} else
				res->comment[index[n]] = NULL;
		}
	}
	res->root = index[tree->root];
	res->size = tree->size;
	res->maxTime = tree->maxTime;
	res->minTime = tree->minTime;
	res->minTimeInt = tree->minTimeInt;
	res->maxTimeInt = tree->maxTimeInt;
	return res;
}

void fillTipBounds(int n, TypeTree *tree, int *start, int *end) {
	if(tree->node[n].child == NOSUCH) {
		start[n] = n;
		end[n] = n;
	} else {
		int c;
		start[n] = INT_MAX;
		end[n] = INT_MIN;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling) {
			fillTipBounds(c, tree, start, end);
			if(start[c]<start[n])
				start[n] = start[c];
			if(end[c]>end[n])
				end[n] = end[c];
		}
	}
}

void fillTimeCumul(int n, TypeTree *tree, double *tcumul) {
	if(tree->node[n].child != NOSUCH) {
		int c;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling) {
			tcumul[c] = tcumul[n]+tree->time[c];
			fillTimeCumul(c, tree, tcumul);
		}
	}
}

void fillCountTips(int n, TypeTree *tree, int *ntips) {
	if(tree->node[n].child == NOSUCH)
		ntips[n] = 1;
	else {
		int c;
		ntips[n] = 0;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling) {
			fillCountTips(c, tree, ntips);
			ntips[n] += ntips[c];
		}
	}
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
	int *parent, n, c, tmp, *path, size, i, *toSetBis, *index, next, out1=NOSUCH, out2=NOSUCH, dec = 0;
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
		out1 = path[size-1];
		if(size == 2) {
			out2 = res->node[path[size-2]].child;
			res->time[path[size-2]] += res->time[res->node[path[size-2]].child];
			res->node[path[size-2]].child = res->node[res->node[path[size-2]].child].child;
		} else {
			out2 = path[size-2];
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
	int n, c, tmp, *parent, *path, size, i;
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
