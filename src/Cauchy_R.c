/*
an object of class "phylo" with the following components:
edge a two-column matrix of mode numeric where each row represents an edge of the
tree; the nodes and the tips are symbolized with numbers; the tips are numbered
1, 2, . . . , and the nodes are numbered after the tips. For each row, the first
column gives the ancestor.
edge.length (optional) a numeric vector giving the lengths of the branches given by edge.
tip.label a vector of mode character giving the names of the tips; the order of the names
in this vector corresponds to the (positive) number in edge.
Nnode the number of (internal) nodes.
node.label (optional) a vector of mode character giving the names of the nodes.
https://informatique-mia.inra.fr/r4ciam/appelC.html
* http://adv-r.had.co.nz/C-interface.html#c-data-structures
* https://cran.r-project.org/doc/manuals/R-ints.html#The-_0027data_0027
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <time.h>
#include "Utils.h"
#include "Tree.h"
#include "SimulateTrait.h"
#include "Cauchy.h"

SEXP getListElement(SEXP list, const char *str) {
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < length(list); i++)
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
        }
    return elmt;
}
	
void inspect(SEXP x) {
  SEXP inspectCall = PROTECT(Rf_lang2(Rf_install("inspect"), x));
  SEXP dotInternalCall = PROTECT(Rf_lang2(Rf_install(".Internal"), inspectCall));
  Rf_eval(dotInternalCall, R_GlobalEnv);
  UNPROTECT(2);
}

TypeTree *Phylo2Tree(SEXP phy) {
	TypeTree *tree;
	int n, i, *parent;
	double *bl;
	SEXP namesTip, edge;
	SEXP dim;
	SEXP rootEdge;
	int nrow, ncol, *tab;

	n = INTEGER(getListElement(phy, "Nnode"))[0];
	bl = REAL(getListElement(phy, "edge.length"));
	edge = getListElement(phy, "edge");
	dim = getAttrib(edge, R_DimSymbol) ;
	rootEdge = getListElement(phy, "root.edge");
	nrow = INTEGER(dim)[0];
	ncol = INTEGER(dim)[1];
	tab = INTEGER(edge);
	tree = (TypeTree *) malloc(sizeof(TypeTree));
	tree->info = NULL;
	tree->comment = NULL;
	tree->parent = NULL;
	tree->size = 0;
	for(i=0; i<nrow*ncol; i++)
		if(tab[i]>tree->size)
			tree->size = tab[i];
	tree->sizeBuf = tree->size;
	tree->node = (TypeNode*) malloc(tree->sizeBuf*sizeof(TypeNode));
	tree->time = (double*) malloc(tree->sizeBuf*sizeof(double));
	tree->name = (char**) malloc(tree->sizeBuf*sizeof(char*));
	for(n=0; n<tree->size; n++) {
		tree->node[n].child = NOSUCH;
		tree->node[n].sibling = NOSUCH;
		tree->name[n] = NULL;
	}
	for(i=0; i<nrow; i++) {
		tree->node[tab[i+nrow]-1].sibling = tree->node[tab[i]-1].child;
		tree->node[tab[i]-1].child = tab[i+nrow]-1;
		tree->time[tab[i+nrow]-1] = bl[i];
	}
	namesTip = getListElement(phy, "tip.label");
	for(i=0; i<LENGTH(namesTip); i++) {
		tree->name[i] = strdpl((char*) CHAR(STRING_ELT(namesTip, i)));
	}
	parent = (int*) malloc(tree->sizeBuf*sizeof(int));
	for(n=0; n<tree->sizeBuf; n++)
		parent[n] = NOSUCH;
	for(n=0; n<tree->size; n++) {
		int c;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling)
			parent[c] = n;
	}
	for(n=0; n<tree->size && parent[n] != NOSUCH; n++)
		;
	free((void*)parent);
	tree->root = n;
	if (isNull(rootEdge)) {
	  tree->time[tree->root] = NO_TIME;
	} else {
          double *rl;
	  rl = REAL(rootEdge);
	  tree->time[tree->root] = rl[0];
	}
	return tree;
}

SEXP printRTree(SEXP phy) {
	FILE *fout;
	if((fout = tmpfile())) {
		TypeTree *tree;
		long numbytes;
		char *buffer;
		tree = Phylo2Tree(phy);
		if (tree->time[tree->root] == NO_TIME)
			tree->time[tree->root] = 0; 
		fprintTreeNewick(fout, tree);
		freeTree(tree);
		fseek(fout, 0L, SEEK_END);
		numbytes = ftell(fout);
		fseek(fout, 0L, SEEK_SET);
		buffer = (char*)calloc(numbytes+1, sizeof(char));
		int size_read = fread(buffer, sizeof(char), numbytes, fout);
		fclose(fout);
		if (size_read != numbytes) 
			error("Temporary file reading failed.");
		buffer[numbytes] = '\0';
		SEXP res = Rf_mkString(buffer);
		free(buffer);
		return res;
	}
	return R_NilValue;
} 

SEXP SimulateTipsCauchy(SEXP treeR, SEXP startR, SEXP dispR) {
	TypeTree *tree;
	double start, disp, *trait;
	TypeCauchyStartSimulData data;
	int n, nTips;
	tree = Phylo2Tree(treeR);
	start = asReal(startR);
	disp = asReal(dispR);
	trait = (double*) malloc(tree->size*sizeof(double));
	data.disp = disp;
	data.start = start;
	GetRNGstate();
	simulTrait(trait, tree, cauchyInitStart, cauchyTransStart, (void*) &data);
	PutRNGstate();
	nTips = 0;
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child == NOSUCH)
			nTips++;
	SEXP res = PROTECT(allocVector(REALSXP, nTips));
	SEXP res_names = PROTECT(res_names = allocVector(STRSXP, nTips));    
	double *rmat = REAL(res);
	nTips = 0;
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child == NOSUCH) {
			rmat[nTips] = trait[n];
			SET_STRING_ELT(res_names, nTips, mkChar(tree->name[n]));
			nTips++;
		}
	setAttrib(res, R_NamesSymbol, res_names);
	UNPROTECT(2);
	freeTree(tree);
	free((void*)trait);
	return res;
}

SEXP getLogDensityTipsCauchy(SEXP treeR, SEXP tipTraitR, SEXP tipNamesR, SEXP startR, SEXP dispR, SEXP typeR, SEXP rootTipR) {
	TypeTree *tree;
	double *trait, dens;
	TypeCauchyInfo *cinf;
	int n, type;
	TypeLexiTree *dict;

	tree = Phylo2Tree(treeR);
	type = asInteger(typeR);
	trait = (double*) malloc(tree->size*sizeof(double));
	for(n=0; n<tree->size; n++)
		trait[n] = DBL_MAX;
	dict = newLexiTree();
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child == NOSUCH && tree->name && tree->name[n]) {
			if(addWordLexi(tree->name[n], n, dict)>=0)
				warning("Warning! duplicate identifier '%s'\n", tree->name[n]);
		}
	for(n=0; n<LENGTH(tipTraitR); n++) {
		int tmp;
		tmp = findWordLexi((char*)CHAR(STRING_ELT(tipNamesR, n)), dict);
		if(tmp != NOSUCH)
			trait[tmp] = REAL(tipTraitR)[n];
	}
	freeLexiTree(dict);
	tree->info = (void*) trait;
	switch(type) {
		case 0:
			cinf = (TypeCauchyInfo*) malloc(tree->size*sizeof(TypeCauchyInfo));
			fillCauchyInfo(tree->root, tree, asReal(dispR), cinf);
			dens = getCauchyLogDensityStem(cinf[tree->root], trait, asReal(startR));
			freeCauchyInfo(tree->root, tree, cinf);
			free((void*)cinf);
			break;
		case 1:
			cinf = (TypeCauchyInfo*) malloc(tree->size*sizeof(TypeCauchyInfo));
			fillCauchyInfo(tree->root, tree, asReal(dispR), cinf);
			dens = getCauchyLogDensityNoStem(cinf[tree->node[tree->root].child], cinf[tree->node[tree->node[tree->root].child].sibling], trait, asReal(startR));
			freeCauchyInfo(tree->root, tree, cinf);
			free((void*)cinf);
			break;
		case 2:
		default:
		{
			TypeTree *rerootedTree;
			  int rootTip;
			  rootTip = asInteger(rootTipR);
			  if(rootTip<0 || rootTip>=tree->size || tree->node[rootTip].child != NOSUCH)
			    for(rootTip=0; rootTip<tree->size && tree->node[rootTip].child != NOSUCH; rootTip++)
			      ;
			rerootedTree = rerootTreeREML(rootTip, tree);
			rerootedTree->info = tree->info;
			cinf = (TypeCauchyInfo*) malloc(tree->size*sizeof(TypeCauchyInfo));
			fillCauchyInfo(rerootedTree->root, rerootedTree, asReal(dispR), cinf);
			dens = getCauchyLogDensityStem(cinf[rerootedTree->root], (double*) rerootedTree->info, ((double*) rerootedTree->info)[rerootedTree->root]);
			freeCauchyInfo(rerootedTree->root, rerootedTree, cinf);
			free((void*)cinf);
			rerootedTree->info = NULL;
			freeTree(rerootedTree);
		}
	}
	free(tree->info);
	tree->info = NULL;
	freeTree(tree);
	return ScalarReal(dens);
}


SEXP getPosteriorLogDensityAncestralCauchy(SEXP nodeR, SEXP tabValR, SEXP treeR, SEXP tipTraitR, SEXP tipNamesR, SEXP startR, SEXP dispR, SEXP typeR) {
	TypeTree *tree;
	double *trait, *dens, *tabVal;
	int n, node, nVal, i, type;
	TypeLexiTree *dict;
	
	type = asInteger(typeR);
	tree = Phylo2Tree(treeR);
	tabVal = REAL(tabValR);
	nVal = LENGTH(tabValR);
	node = asInteger(nodeR);
	trait = (double*) malloc(tree->size*sizeof(double));
	for(n=0; n<tree->size; n++)
		trait[n] = DBL_MAX;
	dict = newLexiTree();
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child == NOSUCH && tree->name && tree->name[n]) {
			if(addWordLexi(tree->name[n], n, dict)>=0)
				warning("Warning! duplicate identifier '%s'\n", tree->name[n]);
		}
	for(n=0; n<LENGTH(tipTraitR); n++) {
		int tmp = findWordLexi((char*)CHAR(STRING_ELT(tipNamesR, n)), dict);
		if(tmp != NOSUCH)
			trait[tmp] = REAL(tipTraitR)[n];
	}
	freeLexiTree(dict);
	tree->info = (void*) trait;
	dens = (double*) malloc(nVal*sizeof(double));
	switch(type) {
		case 0:
			fillCauchyAncestralPosteriorLogDensityStem(node, dens, tabVal, nVal, tree, asReal(dispR), asReal(startR));
			break;
		case 1:
			fillCauchyAncestralPosteriorLogDensityNoStem(node, dens, tabVal, nVal, tree, asReal(dispR), asReal(startR));
			break;
		case 2:
		default:
			fillCauchyAncestralPosteriorLogDensityREML(node, dens, tabVal, nVal, tree, asReal(dispR));
	}
	free(tree->info);
	tree->info = NULL;
	freeTree(tree);

	SEXP densR = PROTECT(allocVector(REALSXP, nVal));
	for(i=0; i<nVal; i++)
		REAL(densR)[i] = dens[i];
	UNPROTECT(1);
	free((void*)dens);
	return densR;
}

SEXP getPosteriorLogDensityIncrementCauchy(SEXP nodeR, SEXP tabValR, SEXP treeR, SEXP tipTraitR, SEXP tipNamesR, SEXP startR, SEXP dispR, SEXP typeR) {
	TypeTree *tree;
	double *trait, *dens, *tabVal;
	int n, node, nVal, i, type;
	TypeLexiTree *dict;
	
	type = asInteger(typeR);
	tree = Phylo2Tree(treeR);
	tabVal = REAL(tabValR);
	nVal = LENGTH(tabValR);
	node = asInteger(nodeR);
	trait = (double*) malloc(tree->size*sizeof(double));
	for(n=0; n<tree->size; n++)
		trait[n] = DBL_MAX;
	dict = newLexiTree();
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child == NOSUCH && tree->name && tree->name[n]) {
			if(addWordLexi(tree->name[n], n, dict)>=0)
				warning("Warning! duplicate identifier '%s'\n", tree->name[n]);
		}
	for(n=0; n<LENGTH(tipTraitR); n++) {
		int tmp = findWordLexi((char*)CHAR(STRING_ELT(tipNamesR, n)), dict);
		if(tmp != NOSUCH)
			trait[tmp] = REAL(tipTraitR)[n];
	}
	freeLexiTree(dict);
	tree->info = (void*) trait;
	dens = (double*) malloc(nVal*sizeof(double));
	switch(type) {
		case 0:
			fillCauchyIncrementPosteriorLogDensityStem(node, dens, tabVal, nVal, tree, asReal(dispR), asReal(startR));
			break;
		case 1:
			fillCauchyIncrementPosteriorLogDensityNoStem(node, dens, tabVal, nVal, tree, asReal(dispR), asReal(startR));
			break;
		case 2:
		default:
			fillCauchyIncrementPosteriorLogDensityREML(node, dens, tabVal, nVal, tree, asReal(dispR));
	}
	free(tree->info);
	tree->info = NULL;
	freeTree(tree);

	SEXP densR = PROTECT(allocVector(REALSXP, nVal));
	for(i=0; i<nVal; i++)
		REAL(densR)[i] = dens[i];
	UNPROTECT(1);
	free((void*)dens);
	return densR;
}
