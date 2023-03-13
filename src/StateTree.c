#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Utils.h"

#include "StateTree.h"

#define INC_STATE_ITEM 50
#define MAX_NAME_SIZE 500

void fixTreeTime(TypeTree *tree) {
	int n;
	if(tree->time == NULL) {
		tree->time = (double*) malloc(tree->sizeBuf*sizeof(double));
		for(n=0; n<tree->sizeBuf; n++)
			tree->time[n] = 1.;
	} else {
		for(n=0; n<tree->size && tree->time[n] != NO_TIME; n++)
			if(tree->time[n] != NO_TIME)
				tree->time[n] = 1.;
	}
	if(tree->minTime == NO_TIME)
		tree->minTime = 0.;	
}
			
		

/*get the status from comments*/
double getState(char *str) {
	int i, ind = 0;
	char tmp[MAX_NAME_SIZE];

	if(str == NULL || strlen(str) < 4)
		return UNKNOWN;
	for(i=3; str[i] != '\0' && (str[i-3]!=':' || str[i-2]!='S' || str[i-1]!='T' || str[i]!='='); i++);
	if(str[i]=='=')
		i++;
	for(; str[i] != '\0' && issep(str[i]); i++);
	for(; str[i] != '\0' && str[i] != ':' && !issep(str[i]); i++) 
		tmp[ind++] = str[i];
	tmp[ind++] = '\0';
	if(ind>1) {
		return atof(tmp);
	} else
		return UNKNOWN;
}

int testLeaves(TypeTree *tree) {
	int n;
	for(n=0; n<tree->size; n++)
		if(isUnknown(n, tree) && tree->node[n].child == -1) {
			Rprintf("Leaf index %d with unknown state\n", n);
//			fprintIdentTimeCommentState(stdout, n, tree, display_both);
			Rprintf("\n");
			return 0;
		}
	return 1;
}

/*return true if the node n is unknown*/
int isUnknown(int n, TypeTree *tree) {
	return ((double*)tree->info)[n] == UNKNOWN;
}

/*from tree to state tree*/
void fillState(TypeTree *tree) {
	int i;
	if(tree->info == NULL) {
		tree->info = (void*) malloc(tree->size*sizeof(double));
		for(i=0; i<tree->size; i++)
			((double*)tree->info)[i] = UNKNOWN;
	}
	for(i=0; i<tree->size; i++)
		((double*)tree->info)[i] = getState(tree->comment[i]);
}

void readStateItem(FILE *f, TypeStateItem **item, int *nitem) {
	int sizeBuf;
	char c;
	
	sizeBuf = INC_STATE_ITEM;
	*item = (TypeStateItem*) malloc(sizeBuf*sizeof(TypeStateItem));
	*nitem = 0;
	do {
		char *tmp;
		int i;
		tmp = (char*) malloc((MAX_NAME_SIZE+1)*sizeof(char));
		skipSeparator(f);
		c = fgetc(f);
		if(c == '\'' || c == '"') {
			c = fgetc(f);
			for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
			if(c == '\'' || c == '"')
				c = fgetc(f);
		} else {
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !IsSeparator(c); i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
		}
		if(i == MAX_NAME_SIZE)
			error("Name too much long...");;
		tmp[i++] = '\0';
		if(i>1) {
			char bof[MAX_NAME_SIZE+1];
			if(*nitem >= sizeBuf) {
				sizeBuf += INC_STATE_ITEM;
				*item = (TypeStateItem*) realloc((void*) *item, sizeBuf*sizeof(TypeStateItem));
			}
			(*item)[*nitem].name = (char *) realloc((void*) tmp, i*sizeof(char));
			skipSeparator(f);
			c = fgetc(f);
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !IsSeparator(c); i++) {
				bof[i] = c;
				c = fgetc(f);
			}
			if(i == MAX_NAME_SIZE)
				error("Number too much long...");;
			bof[i++] = '\0';
			(*item)[*nitem].state = atof(bof);
			(*nitem)++;
		} else {
			free((void*) tmp);
		}
	} while(c != EOF);
	*item = (TypeStateItem*) realloc((void*) *item, (*nitem)*sizeof(TypeStateItem));
}

void freeNameStateList(TypeNameStateList *list) {
	int i;
	for(i=0; i<list->size; i++)
		free((void*)list->name[i]);
	free((void*)list->name);
	free((void*)list->state);
	free((void*)list);
}

void freeIndexStateList(TypeIndexStateList *list) {
	free((void*)list->index);
	free((void*)list->state);
	free((void*)list);
}

TypeNameStateList *readStateList(FILE *f) {
	int sizeBuf;
	char c;
	TypeNameStateList *res;
	
	res = (TypeNameStateList *) malloc(sizeof(TypeNameStateList));
	sizeBuf = INC_STATE_ITEM;
	res->name = (char**) malloc(sizeBuf*sizeof(char*));
	res->state = (double*) malloc(sizeBuf*sizeof(double));
	res->size = 0;
	do {
		char *tmp;
		int i;
		tmp = (char*) malloc((MAX_NAME_SIZE+1)*sizeof(char));
		skipSeparator(f);
		c = fgetc(f);
		if(c == '\'' || c == '"') {
			c = fgetc(f);
			for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
			if(c == '\'' || c == '"')
				c = fgetc(f);
		} else {
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !IsSeparator(c); i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
		}
		if(i == MAX_NAME_SIZE)
			error("Name too much long...");;
		tmp[i++] = '\0';
		removeSpaces(tmp);
		if(strlen(tmp)>0) {
			char bof[MAX_NAME_SIZE+1];
			if(res->size >= sizeBuf) {
				sizeBuf += INC_STATE_ITEM;
				res->name = (char**) realloc((void*) res->name, sizeBuf*sizeof(char*));
				res->state = (double*) realloc((void*) res->state, sizeBuf*sizeof(double));
			}
			res->name[res->size] = (char *) realloc((void*) tmp, i*sizeof(char));
			skipSeparator(f);
			c = fgetc(f);
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !IsSeparator(c); i++) {
				bof[i] = c;
				c = fgetc(f);
			}
			if(i == MAX_NAME_SIZE)
				error("Number too much long...");;
			bof[i++] = '\0';
			res->state[res->size] = atof(bof);
			(res->size)++;
		} else {
			free((void*) tmp);
		}
	} while(c != EOF);
	res->name = (char**) realloc((void*) res->name, res->size*sizeof(char*));
	res->state = (double*) realloc((void*) res->state, res->size*sizeof(double));
	return res;
}

TypeIndexStateList *name2indexStateList(TypeNameStateList *list, TypeTree *tree) {
	TypeLexiTree *dict;
	TypeIndexStateList *res;
	int i;
	res = (TypeIndexStateList *) malloc(sizeof(TypeIndexStateList));
	res->size = 0;
	res->index = (int*) malloc(list->size*sizeof(int));
	res->state = (double*) malloc(list->size*sizeof(double));
	dict = getDictNameTab(tree->name, tree->size);
	for(i=0; i<list->size; i++) {
		int n;
		if((n = findWordLexi(list->name[i], dict))>=0) {
			res->index[res->size] = n;
			res->state[res->size] = list->state[i];
			res->size++;
		}
	}
	res->index = (int*) realloc((void*)res->index, res->size*sizeof(int));
	res->state = (double*) realloc((void*)res->state, res->size*sizeof(double));
	freeLexiTree(dict);
	return res;
}
	
void setState(FILE *f, TypeTree *tree) {
	TypeLexiTree *dict;
	TypeNameStateList *list;
	int i;
	dict = getDictNameTab(tree->name, tree->size);
	list = readStateList(f);
	if(tree->info == NULL) {
		tree->info = (void*) malloc(tree->size*sizeof(double));
		for(i=0; i<tree->size; i++)
			((double*)tree->info)[i] = UNKNOWN;
	}
	for(i=0; i<list->size; i++) {
		int n;
		if((n = findWordLexi(list->name[i], dict))>=0)
			((double*)tree->info)[n] = list->state[i];
	}
	freeNameStateList(list);
	freeLexiTree(dict);
}

/*print ident, time, comment and state of node n*/	
void fprintIdentTimeCommentState(FILE *f, int n, TypeTree *tree, TypeDisplayName display) {
	switch(display) {
		case display_none:
		case display_time_none:
			break;
		case display_name:
		case display_time_name:
			if(tree->name[n] != NULL)
				fprintf(f, "'%s'", tree->name[n]);
			break;
		case display_index:
		case display_time_index:
			fprintf(f, "'%d'", n);
			break;
		case display_both:
		case display_time_both:
		default:
			if(tree->name[n] != NULL)
				fprintf(f, "'%s_", tree->name[n]);
			else
				fprintf(f, "'");
			fprintf(f, "%d'", n);
			break;
	}
	if(display>=display_time_none)
		fprintf(f, ":%lf", tree->time[n]);
	if(tree->comment[n] != NULL)
		fprintf(f, "[%s, ", tree->comment[n]);
	else
		fprintf(f, "[");
	if(isUnknown(n, tree))
		fprintf(f, "ST=?]");
	else
		fprintf(f, "ST=%lf]", ((double*)tree->info)[n]);
}

/*print tab of the states of all nodes*/	
void fprintStateTreeList(FILE *f, TypeTree *tree, TypeDisplayName display) {
	int n;
	for(n=0; n<tree->size; n++) {
		fprintf(f, "%d\t", n);
		fprintIdentTimeCommentState(f, n, tree, display);
		fprintf(f, "\t%.1lf\n", ((double*)tree->info)[n]);
	}
}

/*print tab of the states of all nodes*/	
void fprintStateList(FILE *f, TypeNameStateList *list) {
	int n;
	for(n=0; n<list->size; n++) {
		fprintf(f, "'%s'\t%lE\n", list->name[n], list->state[n]);
	}
}
/*print tab of the states of leaves*/	
void fprintStateLeavesList(FILE *f, TypeTree *tree, TypeDisplayName display) {
	int n;
	for(n=0; n<tree->size; n++) {
		if(tree->node[n].child<0) {
			if(tree->name[n] != NULL)
				fprintf(f, "%s", tree->name[n]);
			fprintf(f, "\t%.1lf\n", ((double*)tree->info)[n]);
		}
	}
}
/*print tree in newick format*/	
void fprintStateTree(FILE *f, TypeTree *tree, TypeDisplayName display) {
	if(tree->size<=0)
		return;
	if(tree->node[tree->root].child >= 0) {
		int tmp = tree->node[tree->root].child;
		fprintf(f, "(");
		fprintStateNode(f, tmp, tree, display);
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
			fprintf(f, ", ");
			fprintStateNode(f, tmp, tree, display);
		}
		fprintf(f, ")");
	}
	fprintIdentTimeCommentState(f, tree->root, tree, display);
	fprintf(f, ";\n");
}

/*print node in newick format*/	
void fprintStateNode(FILE *f, int n, TypeTree *tree, TypeDisplayName display) {
	if(tree->node[n].child >= 0) {
		int tmp = tree->node[n].child;
		fprintf(f, "(");
		fprintStateNode(f, tmp, tree, display);
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
			fprintf(f, ", ");
			fprintStateNode(f, tmp, tree, display);
		}
		fprintf(f, ")");
	}
	fprintIdentTimeCommentState(f, n, tree, display);
}

