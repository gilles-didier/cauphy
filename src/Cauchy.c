#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include "Utils.h"
#include "Kahan.h"
#include "Cauchy.h"

#define Cauchy_MIN(x,y) ((x)<(y)?(x):(y))

static int compare_child(double ta, double va, double tb, double vb);

typedef struct CHILD_SINGLE_TMP {
	int loc, gen;
} TypeChildSingleTmp;

typedef struct CHILD_IDENT_TMP {
	int locl, locr, gen;
} TypeChildIdentTmp;

int compare_child(double ta, double va, double tb, double vb) {
	if(ta > tb)
		return 1;
	else {
		if(ta < tb)
			return -1;
		else {
			if(va > vb)
				return 1;
			else {
				if(va < vb)
					return -1;
				else 
					return 0;
			}
		}
	}
}

double getCauchyLogDensityStandard(double x,  double d) {
	return log(d)-log(M_PI)-log(x*x+d*d);
}

void freeCauchyInfo(int n, TypeTree *tree, TypeCauchyInfo *cinf) {
	if(tree->node[n].child != NOSUCH) {
		freeCauchyInfo(tree->node[n].child, tree, cinf);
		freeCauchyInfo(tree->node[tree->node[n].child].sibling, tree, cinf);
	}
	free((void*)cinf[n].val);
	free((void*)cinf[n].time);
	free((void*)cinf[n].A);
}
		
void fillCauchyInfo(int n, TypeTree *tree, double param, TypeCauchyInfo *cinf) {
  	if(tree->node[n].child == NOSUCH) {
		if(isnan(((double*)tree->info)[n]) || !isfinite(((double*)tree->info)[n])) {
			if(tree->name && tree->name[n])
				error("Value of %s is not a valid number (%le).\n", tree->name[n], ((double*)tree->info)[n]);
			else
				error("Value of %s is not a valid number (%le).\n", n, ((double*)tree->info)[n]);
		}
		cinf[n].sizeChild = 1;
		cinf[n].time = (double*) malloc(cinf[n].sizeChild*sizeof(double));
		cinf[n].val = (double*) malloc(cinf[n].sizeChild*sizeof(double));
		cinf[n].A = (TypeComplex*) malloc(cinf[n].sizeChild*sizeof(TypeComplex));
		cinf[n].time[0] = param*tree->time[n];
		cinf[n].val[0] = ((double*)tree->info)[n];
		cinf[n].A[0] = getComplex(1., 0.);
	} else {
		int cl, cr, i, j, k, size, size_sr, size_sl, size_id;
		TypeChildSingleTmp *sl, *sr;
		TypeChildIdentTmp *id;
		cl = tree->node[n].child;
		cr = tree->node[tree->node[n].child].sibling;
		if(cr == NOSUCH) {
			if(tree->name && tree->name[n])
				error("Node %s has a single child\n", tree->name[n]);
			else
				error("Node n %d has a single child\n", n);
		}
		fillCauchyInfo(cl, tree, param, cinf);
		fillCauchyInfo(cr, tree, param, cinf);
		sl = (TypeChildSingleTmp*) malloc(cinf[cl].sizeChild*sizeof(TypeChildSingleTmp));
		sr = (TypeChildSingleTmp*) malloc(cinf[cr].sizeChild*sizeof(TypeChildSingleTmp));
		id = (TypeChildIdentTmp*) malloc(Cauchy_MIN(cinf[cl].sizeChild, cinf[cr].sizeChild)*sizeof(TypeChildIdentTmp));
		size = 0; size_sr = 0; size_sl = 0; size_id = 0;
		i = 0; j = 0;
		while(i<cinf[cl].sizeChild || j<cinf[cr].sizeChild) {
			if(i>=cinf[cl].sizeChild) {
				sr[size_sr].loc = j;
				sr[size_sr].gen = size;
				j++; size_sr++; size++;
			} else {
				if(j>=cinf[cr].sizeChild) {
					sl[size_sl].loc = i;
					sl[size_sl].gen = size;
					i++; size_sl++; size++;
				} else {
					switch(compare_child(cinf[cl].time[i], cinf[cl].val[i], cinf[cr].time[j], cinf[cr].val[j])) {
						case -1:
							sl[size_sl].loc = i;
							sl[size_sl].gen = size;
							i++; size_sl++; size++;
							break;
						case 0:
							id[size_id].locl = i;
							id[size_id].locr = j;
							id[size_id].gen = size;
							i++, j++; size_id++; size++;
							break;
						case 1:
							sr[size_sr].loc = j;
							sr[size_sr].gen = size;
							j++; size_sr++; size++;
							break;
						default:
						;
					}
				}
			}
		}
		cinf[n].sizeChild = size;
		cinf[n].time = (double*) malloc(cinf[n].sizeChild*sizeof(double));
		cinf[n].val = (double*) malloc(cinf[n].sizeChild*sizeof(double));
		cinf[n].A = (TypeComplex*) malloc(cinf[n].sizeChild*sizeof(TypeComplex));
		for(i=0; i<size_sl; i++) {
			cinf[n].val[sl[i].gen] = cinf[cl].val[sl[i].loc];
			cinf[n].time[sl[i].gen] = cinf[cl].time[sl[i].loc]+param*tree->time[n];
			cinf[n].A[sl[i].gen] = getComplex(0., 0.);
			TypeComplexKahan ka;
			initComplexKahan(&ka);
			for(j=0; j<cinf[cr].sizeChild; j++) {
				double x;
				x = cinf[cr].val[j]-cinf[n].val[sl[i].gen];
				sumComplexKahan(divComplex(conjComplex(cinf[cr].A[j]), getComplex(cinf[cr].time[j]+cinf[cl].time[sl[i].loc], -x)), &ka);
				sumComplexKahan(divComplex(cinf[cr].A[j], getComplex(cinf[cr].time[j]-cinf[cl].time[sl[i].loc], x)), &ka);
			}
			cinf[n].A[sl[i].gen] = prodComplex(totalComplexKahan(&ka), cinf[cl].A[sl[i].loc]);
		}
		for(j=0; j<size_sr; j++) {
			cinf[n].val[sr[j].gen] = cinf[cr].val[sr[j].loc];
			cinf[n].time[sr[j].gen] = cinf[cr].time[sr[j].loc]+param*tree->time[n];
			cinf[n].A[sr[j].gen] = getComplex(0., 0.);
			TypeComplexKahan ka;
			initComplexKahan(&ka);
			for(i=0; i<cinf[cl].sizeChild; i++) {
				double x;
				x = cinf[cl].val[i]-cinf[cr].val[sr[j].loc];
				sumComplexKahan(divComplex(conjComplex(cinf[cl].A[i]), getComplex(cinf[cl].time[i]+cinf[cr].time[sr[j].loc], -x)), &ka);
				sumComplexKahan(divComplex(cinf[cl].A[i], getComplex(cinf[cl].time[i]-cinf[cr].time[sr[j].loc], x)), &ka);
			}
			cinf[n].A[sr[j].gen] = prodComplex(totalComplexKahan(&ka), cinf[cr].A[sr[j].loc]);
		}
		for(k=0; k<size_id; k++) {
			cinf[n].val[id[k].gen] = cinf[cr].val[id[k].locr];
			cinf[n].time[id[k].gen] = cinf[cr].time[id[k].locr]+param*tree->time[n];
			cinf[n].A[id[k].gen] = getComplex(0., 0.);
			TypeComplexKahan ka;
			initComplexKahan(&ka);
			for(j=0; j<id[k].locr; j++) {
				double x;
				x = cinf[cr].val[j]-cinf[cl].val[id[k].locl];
				sumComplexKahan(divComplex(conjComplex(cinf[cr].A[j]), getComplex(cinf[cr].time[j]+cinf[cl].time[id[k].locl], -x)), &ka);
				sumComplexKahan(divComplex(cinf[cr].A[j], getComplex(cinf[cr].time[j]-cinf[cl].time[id[k].locl], x)), &ka);
			}
			for(j=id[k].locr+1; j<cinf[cr].sizeChild; j++) {
				double x;
				x = cinf[cr].val[j]-cinf[cl].val[id[k].locl];
				sumComplexKahan(divComplex(conjComplex(cinf[cr].A[j]), getComplex(cinf[cr].time[j]+cinf[cl].time[id[k].locl], -x)), &ka);
				sumComplexKahan(divComplex(cinf[cr].A[j], getComplex(cinf[cr].time[j]-cinf[cl].time[id[k].locl], x)), &ka);
			}
			cinf[n].A[id[k].gen] = prodComplex(totalComplexKahan(&ka), cinf[cl].A[id[k].locl]);
			initComplexKahan(&ka);
			for(i=0; i<id[k].locl; i++) {
				double x;
				x = cinf[cl].val[i]-cinf[cr].val[id[k].locr];
				sumComplexKahan(divComplex(conjComplex(cinf[cl].A[i]), getComplex(cinf[cl].time[i]+cinf[cr].time[id[k].locr], -x)), &ka);
				sumComplexKahan(divComplex(cinf[cl].A[i], getComplex(cinf[cl].time[i]-cinf[cr].time[id[k].locr], x)), &ka);
			}
			for(i=id[k].locl+1; i<cinf[cl].sizeChild; i++) {
				double x;
				x = cinf[cl].val[i]-cinf[cr].val[id[k].locr];
				sumComplexKahan(divComplex(conjComplex(cinf[cl].A[i]), getComplex(cinf[cl].time[i]+cinf[cr].time[id[k].locr], -x)), &ka);
				sumComplexKahan(divComplex(cinf[cl].A[i], getComplex(cinf[cl].time[i]-cinf[cr].time[id[k].locr], x)), &ka);
			}
			cinf[n].A[id[k].gen] = addComplex(cinf[n].A[id[k].gen], prodComplex(totalComplexKahan(&ka), cinf[cr].A[id[k].locr]));
			cinf[n].A[id[k].gen] = addComplex(cinf[n].A[id[k].gen], getComplex((getRealComplex(cinf[cl].A[id[k].locl])*getRealComplex(cinf[cr].A[id[k].locr])-getImagComplex(cinf[cl].A[id[k].locl])*getImagComplex(cinf[cr].A[id[k].locr]))/cinf[cl].time[id[k].locl], 0.));
		}
		free((void*) sl);
		free((void*) sr);
		free((void*) id);
	}
}

double getCauchyLogDensityStem(TypeCauchyInfo cinf, double start) {
	int i, *signA;
	double res, *allLogA;
	signA = (int*) malloc(cinf.sizeChild*sizeof(int));
	allLogA = (double*) malloc(cinf.sizeChild*sizeof(double));
	TypeComplexKahan ka;
	initComplexKahan(&ka);
	for(i=0; i<cinf.sizeChild; i++) {
		double costheta;
		TypeComplex logd1;
		logd1 = getComplex(cinf.time[i], cinf.val[i]-start);
		costheta = cos(getImagLogComplex(cinf.A[i])-getImagLogComplex(logd1));
		signA[i] = sign(costheta);
		allLogA[i] = getRealLogComplex(cinf.A[i])-getRealLogComplex(logd1)+log(fabs(costheta));
	}
	res = log(2.) - ((double)cinf.sizeChild) * log(2 * M_PI)+ sumSignedLogTableKahan(allLogA, signA, cinf.sizeChild);
	free((void*)signA);
	free((void*)allLogA);
	return res;
}

double getCauchyLogDensityNoStem(TypeCauchyInfo cinfL, TypeCauchyInfo cinfR, double start) {
	return getCauchyLogDensityStem(cinfL, start)+getCauchyLogDensityStem(cinfR, start);
}

void fillCauchyAncestralPosteriorLogDensityStem(int node, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp, double start) {
	if(tree->node[node].child == NOSUCH) {
		int j;
		for(j=0; j<nVal; j++)
			if(fabs(tabVal[j] - ((double*)tree->info)[node]) < 0.00001)
				dens[j] = 0.;
			else
				dens[j] = log(0.);
	} else {
		TypeCauchyInfo *cinf;
		int j;
		double densRef;
		cinf = (TypeCauchyInfo*) malloc(tree->sizeBuf*sizeof(TypeCauchyInfo));
		fillCauchyInfo(tree->root, tree, disp, cinf);
		densRef = getCauchyLogDensityStem(cinf[tree->root], start);
		if(node == tree->root) {
			for(j=0; j<nVal; j++)
				dens[j] = getCauchyLogDensityStem(cinf[tree->node[node].child], tabVal[j])
				+getCauchyLogDensityStem(cinf[tree->node[tree->node[node].child].sibling], tabVal[j])
				+getCauchyLogDensityStandard(tabVal[j]-start, tree->time[node]*disp)
				-densRef;		
			freeCauchyInfo(tree->root, tree, cinf);
		} else {
			TypeTree *treeC;
			int childSave;
			childSave = tree->node[node].child;
			tree->node[node].child = NOSUCH;
			treeC = rerootTreeStem(node, tree);
			treeC->info = tree->info;
			tree->node[node].child = childSave;
			((double*)tree->info)[tree->root] = start;
			freeCauchyInfo(tree->root, tree, cinf);
			fillCauchyInfo(tree->node[node].child, tree, disp, cinf);
			fillCauchyInfo(tree->node[tree->node[node].child].sibling, tree, disp, cinf);
			fillCauchyInfo(treeC->root, treeC, disp, cinf);
			for(j=0; j<nVal; j++) {
				dens[j] = getCauchyLogDensityStem(cinf[tree->node[node].child], tabVal[j])
				+getCauchyLogDensityStem(cinf[tree->node[tree->node[node].child].sibling], tabVal[j])
				+getCauchyLogDensityStem(cinf[treeC->root], tabVal[j])
				-densRef;
			}
			freeCauchyInfo(tree->node[node].child, tree, cinf);
			freeCauchyInfo(tree->node[tree->node[node].child].sibling, tree, cinf);
			freeCauchyInfo(treeC->root, treeC, cinf);
			treeC->info = NULL;
			freeTree(treeC);
		}
		free((void*)cinf);
	}
}

void fillCauchyAncestralPosteriorLogDensityNoStem(int node, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp, double start) {
	if(node == tree->root) {
		int j;
		for(j=0; j<nVal; j++)
			if(fabs(tabVal[j] - start) < 0.00001)
				dens[j] = 0.;
			else
				dens[j] = log(0.);
	} else {
		int rootSave;
		rootSave = tree->root;
		tree->root = findSide(node, tree);
		fillCauchyAncestralPosteriorLogDensityStem(node, dens, tabVal, nVal, tree, disp, start);
		tree->root = rootSave;
	}
}

void fillCauchyAncestralPosteriorLogDensityREML(int node, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp) {
	TypeTree *rerootedTree;
	double start;
	int rootTip;
	if(node == tree->root)
		for(rootTip=node; tree->node[rootTip].child != NOSUCH; rootTip = tree->node[rootTip].child)
			;
	else {
		int *parent;
		parent = getParent(tree);
		for(rootTip=(node!=tree->node[parent[node]].child)?tree->node[parent[node]].child:tree->node[tree->node[parent[node]].child].sibling; tree->node[rootTip].child != NOSUCH; rootTip = tree->node[rootTip].child)
			;
		if(parent[rootTip] == tree->root && parent[node] == tree->root)
			node = rootTip;
		free((void*)parent);
	}	
	rerootedTree = rerootTreeREML(rootTip, tree);
	rerootedTree->info = tree->info;
	if(rerootedTree->size > 1)
		start = ((double*) rerootedTree->info)[rerootedTree->root];
	else
		start = (tree->node[tree->root].child == rootTip)?((double*) rerootedTree->info)[tree->node[tree->node[tree->root].child].sibling]:((double*) rerootedTree->info)[tree->node[tree->root].child];
	if(node != tree->root)
		fillCauchyAncestralPosteriorLogDensityStem(node, dens, tabVal, nVal, rerootedTree, disp, start);
	else {
		int i;
		double densRef;
		TypeCauchyInfo *cinf;
		cinf = (TypeCauchyInfo*) malloc(tree->size*sizeof(TypeCauchyInfo));
		fillCauchyInfo(rerootedTree->root, rerootedTree, disp, cinf);
		densRef = getCauchyLogDensityStem(cinf[rerootedTree->root],  start);
		freeCauchyInfo(rerootedTree->root, rerootedTree, cinf);
		fillCauchyInfo(tree->root, tree, disp, cinf);
		for(i=0; i<nVal; i++)
			dens[i] = getCauchyLogDensityNoStem(cinf[tree->node[tree->root].child], cinf[tree->node[tree->node[tree->root].child].sibling], tabVal[i])-densRef;
		freeCauchyInfo(tree->root, tree, cinf);
		free((void*)cinf);		
	}
	rerootedTree->info = NULL;
	freeTree(rerootedTree);
}

void fillCauchyIncrementPosteriorLogDensityStem(int node, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp, double start) {	
	if(tree->node[node].child == NOSUCH || node == tree->root) {
		int j;
		if(tree->node[node].child == NOSUCH && node == tree->root) {
			for(j=0; j<nVal; j++)
				if(fabs(tabVal[j] - (((double*)tree->info)[node]-start)) < 0.00001)
					dens[j] = 0.;
				else
					dens[j] = log(0.);
		} else {
			double densRef;
			TypeCauchyInfo *cinf;
			if(tree->node[node].child == NOSUCH) {
				TypeTree *rerootedTree;
				((double*)tree->info)[tree->root] = start;
				rerootedTree = rerootTreeStem(node, tree);
				rerootedTree->info = tree->info;
				start = ((double*)tree->info)[node];
				cinf = (TypeCauchyInfo*) malloc(rerootedTree->sizeBuf*sizeof(TypeCauchyInfo));
				fillCauchyInfo(rerootedTree->root, rerootedTree, disp, cinf);
				densRef = getCauchyLogDensityStem(cinf[rerootedTree->root], start);
				for(j=0; j<nVal; j++)
					dens[j] = getCauchyLogDensityStem(cinf[rerootedTree->node[rerootedTree->root].child], start-tabVal[j]) +
						getCauchyLogDensityStem(cinf[rerootedTree->node[rerootedTree->node[rerootedTree->root].child].sibling], start-tabVal[j])+
						getCauchyLogDensityStandard(tabVal[j], disp*tree->time[node]) - densRef;
				freeCauchyInfo(rerootedTree->root, rerootedTree, cinf);
				free((void*)cinf);
				rerootedTree->info = NULL;
				freeTree(rerootedTree);
			} else {
				cinf = (TypeCauchyInfo*) malloc(tree->sizeBuf*sizeof(TypeCauchyInfo));
				fillCauchyInfo(tree->root, tree, disp, cinf);
				densRef = getCauchyLogDensityStem(cinf[tree->root], start);
				for(j=0; j<nVal; j++)
					dens[j] = getCauchyLogDensityStem(cinf[tree->node[tree->root].child], start+tabVal[j]) +
						getCauchyLogDensityStem(cinf[tree->node[tree->node[tree->root].child].sibling], start+tabVal[j])+
						getCauchyLogDensityStandard(tabVal[j], disp*tree->time[node]) - densRef;
				freeCauchyInfo(tree->root, tree, cinf);
				free((void*)cinf);
			}
		}
	} else {
		int i, j, *tips, nTips;
		double densRef, timeSave;
		void *infoSave;
		TypeCauchyInfo *cinf;
		tips = (int*) malloc(tree->size*sizeof(int));
		nTips = fillTips(node, tree, tips);
		cinf = (TypeCauchyInfo*) malloc(tree->sizeBuf*sizeof(TypeCauchyInfo));
		fillCauchyInfo(tree->root, tree, disp, cinf);
		densRef = getCauchyLogDensityStem(cinf[tree->root], start);
		freeCauchyInfo(tree->root, tree, cinf);
		infoSave = tree->info;
		tree->info = (double*) malloc(tree->size*sizeof(double));
		for(i=0; i<tree->size; i++)
			((double*)tree->info)[i] = ((double*)infoSave)[i];
		timeSave = tree->time[node];
		tree->time[node] = 0.;
		for(j=0; j<nVal; j++) {
			for(i=0; i<nTips; i++)
				((double*)tree->info)[tips[i]] = ((double*)infoSave)[tips[i]]-tabVal[j];
			fillCauchyInfo(tree->root, tree, disp, cinf);
			dens[j] = getCauchyLogDensityStem(cinf[tree->root], start) + getCauchyLogDensityStandard(tabVal[j], disp*timeSave) - densRef;
			freeCauchyInfo(tree->root, tree, cinf);
		}
		free((void*)tips);
		free((void*)tree->info);
		tree->info = infoSave;
		tree->time[node] = timeSave;
		free((void*)cinf);
	}
}

void fillCauchyIncrementPosteriorLogDensityNoStem(int node, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp, double start) {
	if(node != tree->root) {
		int rootSave;
		rootSave = tree->root;
		tree->root = findSide(node, tree);
		fillCauchyIncrementPosteriorLogDensityStem(node, dens, tabVal, nVal, tree, disp, start);
		tree->root = rootSave;
	} else
		error("No stem trees do not have a branch ending with the root\n");
}

void fillCauchyIncrementPosteriorLogDensityREML(int node, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp) {
	if(node != tree->root) {
		if(node != tree->node[tree->root].child && node != tree->node[tree->node[tree->root].child].sibling) {
			TypeTree *rerootedTree;
			int rootTip, *parent;
			parent = getParent(tree);
			for(rootTip=(node!=tree->node[parent[node]].child)?tree->node[parent[node]].child:tree->node[tree->node[parent[node]].child].sibling; tree->node[rootTip].child != NOSUCH; rootTip = tree->node[rootTip].child)
				;
			free((void*)parent);
			rerootedTree = rerootTreeREML(rootTip, tree);
			rerootedTree->info = tree->info;
			fillCauchyIncrementPosteriorLogDensityStem(node, dens, tabVal, nVal, rerootedTree, disp, ((double*) rerootedTree->info)[rerootedTree->root]);
			rerootedTree->info = NULL;
			freeTree(rerootedTree);
		} else {
			TypeCauchyInfo *cinf;
			TypeTree *rerootedTree;
			double densRef;
			int rootTip, other;
			cinf = (TypeCauchyInfo*) malloc(tree->size*sizeof(TypeCauchyInfo));
			other = (node == tree->node[tree->root].child)? tree->node[tree->node[tree->root].child].sibling : tree->node[tree->root].child;
			for(rootTip=other; tree->node[rootTip].child != NOSUCH; rootTip = tree->node[rootTip].child)
				;
			rerootedTree = rerootTreeREML(rootTip, tree);
			rerootedTree->info = tree->info;
			fillCauchyInfo(rerootedTree->root, rerootedTree, disp, cinf);
			if(rerootedTree->size > 1)
				densRef = getCauchyLogDensityStem(cinf[rerootedTree->root], ((double*)rerootedTree->info)[rerootedTree->root]);
			else
				densRef = getCauchyLogDensityStem(cinf[rerootedTree->root], ((double*)tree->info)[node]);				
			freeCauchyInfo(rerootedTree->root, rerootedTree, cinf);
			if(tree->node[node].child == NOSUCH) {
				int j;
				fillCauchyInfo(other, tree, disp, cinf);
				for(j=0; j<nVal; j++)
					dens[j] = getCauchyLogDensityStem(cinf[other], ((double*)tree->info)[node]-tabVal[j]) + getCauchyLogDensityStandard(tabVal[j], disp*tree->time[node]) - densRef;
				freeCauchyInfo(other, tree, cinf);
			} else {
				int i, j, *tips, nTips;
				double time = tree->time[node];
				if(rootTip == other)
					node = other;
				tips = (int*) malloc(tree->size*sizeof(int));
				nTips = fillTips(node, tree, tips);
				rerootedTree->info = malloc(rerootedTree->size*sizeof(double));
				rerootedTree->time[node] = tree->time[other];
				for(j=0; j<rerootedTree->size; j++)
					((double*)rerootedTree->info)[j] = ((double*)tree->info)[j];
				for(j=0; j<nVal; j++) {
					for(i=0; i<nTips; i++)
						((double*)rerootedTree->info)[tips[i]] = ((double*)tree->info)[tips[i]]-tabVal[j];
					fillCauchyInfo(rerootedTree->root, rerootedTree, disp, cinf);
					dens[j] = getCauchyLogDensityStem(cinf[rerootedTree->root], ((double*)rerootedTree->info)[rerootedTree->root]) + getCauchyLogDensityStandard(tabVal[j], disp*time) - densRef;
					freeCauchyInfo(rerootedTree->root, rerootedTree, cinf);
				}
				free((void*)tips);
				free(rerootedTree->info);
			}
			free((void*)cinf);
			rerootedTree->info = NULL;
			freeTree(rerootedTree);
		}		
	} else
		error("Cannot compute the density of the increment of the branch ending with root in the REML case\n");
}
