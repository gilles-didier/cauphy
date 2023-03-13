/*
    'tdalp' computes the time-dependent-asymmetric-linear-parsimony parametric reconstruction of a phylogenetic tree.

    Copyright (C) 2016  Gilles DIDIER

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#ifndef StateTreeF
#define StateTreeF

#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Tree.h"

#define UNKNOWN_STRING "?"



typedef struct TMP_STATE_ITEM {
	char *name;
	double state;
} TypeStateItem;


typedef struct NAME_STATE_LIST {
	char **name;
	double *state;
	int size;
} TypeNameStateList;

typedef struct INDEX_STATE_LIST {
	int size, *index;
	double *state;
} TypeIndexStateList;


void fixTreeTime(TypeTree *tree);
TypeNameStateList *readStateList(FILE *f);
void freeNameStateList(TypeNameStateList *list);
void freeIndexStateList(TypeIndexStateList *list);
void readStateItem(FILE *f, TypeStateItem **item, int *nitem);
TypeNameStateList *readStateList(FILE *f);
double getState(char *str);
void setState(FILE *f, TypeTree *tree);
int isUnknown(int n, TypeTree *tree);
int countLeaves(TypeTree *tree);
void fillState(TypeTree *tree);
void fprintStateTree(FILE *f, TypeTree *tree, TypeDisplayName display);
void fprintStateNode(FILE *f, int n, TypeTree *tree, TypeDisplayName display);
void fprintIdentTimeCommentState(FILE *f, int n, TypeTree *tree, TypeDisplayName display);
/*print tab of the states of all nodes*/	
void fprintStateTreeList(FILE *f, TypeTree *tree, TypeDisplayName display);
/*print tab of the states of all nodes*/	
void fprintStateList(FILE *f, TypeNameStateList *list);
void fprintStateLeavesList(FILE *f, TypeTree *tree, TypeDisplayName display);
int testLeaves(TypeTree *tree);
TypeIndexStateList *name2indexStateList(TypeNameStateList *list, TypeTree *tree);

#endif
