#include <stdlib.h>
#include <string.h>
#include "Utils.h"

#define INC_SIZE_DICT 100
#define MAX_DICT_LENGTH 1000

char *strdpl(char *src) {
    if(src != NULL) {
        char *res;
        res = (char*) malloc((strlen(src)+1)*sizeof(char));
        strcpy(res, src);
        return res;
    } else
        return NULL;
}

int findWordLexi(char *w, TypeLexiTree *dict) {
    int i, n;
    n = dict->root;
    i = 0;
    while(i<MAX_DICT_LENGTH && w[i] != '\0') {
        int c;
        c = dict->node[n].child;
        while(c>=0 && dict->node[c].symbol<w[i])
            c = dict->node[c].sibling;
        if(c>=0 && dict->node[c].symbol == w[i]) {
            n = c;
            i++;
        } else
            return -1;
    }
    if(dict->node[n].child >= 0 && dict->node[dict->node[n].child].symbol == '\0')
        return dict->node[dict->node[n].child].child;
    else
        return -1;
}

/*add word w in dict. If w is already in, then it returns its index (stored in the child field of a leaf - labelled by '\0'), -1 otherwise*/
int addWordLexi(char *w, int index, TypeLexiTree *dict) {
    int i, n;
    n = dict->root;
    i = 0;
    if(w == NULL)
		return -1;
    while(i<MAX_DICT_LENGTH && w[i] != '\0') {
        int *prec, c;
        prec = &(dict->node[n].child);
        c = dict->node[n].child;
        while(c>=0 && dict->node[c].symbol<w[i]) {
            prec = &(dict->node[c].sibling);
            c = dict->node[c].sibling;
        }
        if(c>=0 && dict->node[c].symbol == w[i]) {
            n = c;
            i++;
        } else {
            *prec = dict->size;
            if(dict->size >= dict->sizeBuf) {
                dict->sizeBuf += INC_SIZE_DICT;
                dict->node = (TypeLexiNode*) realloc((void*)dict->node, dict->sizeBuf*sizeof(TypeLexiNode));
            }
            initLexiNode(w[i], &(dict->node[dict->size]));
            dict->node[dict->size].sibling = c;
            n = dict->size;
            dict->size++;
            i++;
            while(i<MAX_DICT_LENGTH && w[i] != '\0') {
                if(dict->size >= dict->sizeBuf) {
                    dict->sizeBuf += INC_SIZE_DICT;
                    dict->node = (TypeLexiNode*) realloc((void*)dict->node, dict->sizeBuf*sizeof(TypeLexiNode));
                }
                initLexiNode(w[i], &(dict->node[dict->size]));
                dict->node[n].child = dict->size;
                n = dict->node[n].child;
                dict->size++;
                i++;
            }
        }
    }
    if(dict->node[n].child >= 0 && dict->node[dict->node[n].child].symbol == '\0')
        return dict->node[dict->node[n].child].child;
    if(dict->size >= dict->sizeBuf) {
        dict->sizeBuf += INC_SIZE_DICT;
        dict->node = (TypeLexiNode*) realloc((void*)dict->node, dict->sizeBuf*sizeof(TypeLexiNode));
    }
    initLexiNode('\0', &(dict->node[dict->size]));
    dict->node[dict->size].sibling = dict->node[n].child;
    dict->node[dict->size].child = index;
    dict->node[n].child = dict->size;
    dict->size++;
    return -1;
}

void initLexiNode(char symbol, TypeLexiNode *n) {
    n->symbol = symbol;
    n->child = -1;
    n->sibling = -1;
}

TypeLexiTree *newLexiTree() {
    TypeLexiTree* dict;
    dict = (TypeLexiTree*) malloc(sizeof(TypeLexiTree));
    dict->sizeBuf = INC_SIZE_DICT;
    dict->node = (TypeLexiNode*) malloc(dict->sizeBuf*sizeof(TypeLexiNode));
    dict->root = 0;
    initLexiNode(0,&(dict->node[dict->root]));
    dict->size = 1;
    return dict;
}

void freeLexiTree(TypeLexiTree *dict) {
    if(dict->node != NULL)
        free((void*)dict->node);
    free((void*) dict);
}











