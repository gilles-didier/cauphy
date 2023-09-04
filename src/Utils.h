#ifndef UtilsF
#define UtilsF

/*if symbol ==-1 then child contains an index*/
typedef struct LEXI_NODE {
    char symbol;
    int child, sibling;
} TypeLexiNode;

typedef struct LEXI_TREE {
    TypeLexiNode *node;
    int root, size, sizeBuf;
} TypeLexiTree;

#ifndef M_PI
#define M_PI acos(-1.)
#endif


#include <R.h>

#ifdef __cplusplus
extern "C" {
#endif
char *strdpl(char *src);
int findWordLexi(char *w, TypeLexiTree *dict);
int addWordLexi(char *w, int index, TypeLexiTree *dict);
void initLexiNode(char symbol, TypeLexiNode *n);
TypeLexiTree *newLexiTree(void);
void freeLexiTree(TypeLexiTree *dict);

#ifdef __cplusplus
}
#endif

#endif
