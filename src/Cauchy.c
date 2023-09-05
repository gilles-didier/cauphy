#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include "Utils.h"
#include "Kahan.h"
#include "Cauchy.h"

double getCauchyLogDensityStandard(double x,  double d) {
	return log(d)-log(M_PI)-log(x*x+d*d);
}

void freeCauchyInfo(int n, TypeTree *tree, TypeCauchyInfo *cinf) {
	if(tree->node[n].child != NOSUCH) {
		freeCauchyInfo(tree->node[n].child, tree, cinf);
		freeCauchyInfo(tree->node[tree->node[n].child].sibling, tree, cinf);
	}
	free((void*)cinf[n].child);
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
		cinf[n].child = (int*) malloc(cinf[n].sizeChild*sizeof(int));
		cinf[n].time = (double*) malloc(cinf[n].sizeChild*sizeof(double));
		cinf[n].A = (TypeComplex*) malloc(cinf[n].sizeChild*sizeof(TypeComplex));
		cinf[n].child[0] = n;
		cinf[n].time[0] = param*tree->time[n];
		cinf[n].A[0] = getComplex(1., 0.);
	} else {
		int cl, cr, i, j;
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
		cinf[n].sizeChild = cinf[cl].sizeChild+cinf[cr].sizeChild;
		cinf[n].child = (int*) malloc(cinf[n].sizeChild*sizeof(int));
		cinf[n].time = (double*) malloc(cinf[n].sizeChild*sizeof(double));
		cinf[n].A = (TypeComplex*) malloc(cinf[n].sizeChild*sizeof(TypeComplex));
		for(i=0; i<cinf[cl].sizeChild; i++) {
			cinf[n].child[i] = cinf[cl].child[i];
			cinf[n].time[i] = cinf[cl].time[i]+param*tree->time[n];
			cinf[n].A[i] = getComplex(0., 0.);
			TypeComplexKahan ka;
			initComplexKahan(&ka);
			for(j=0; j<cinf[cr].sizeChild; j++) {
				double x;
				x = ((double*)tree->info)[cinf[cr].child[j]]-((double*)tree->info)[cinf[cl].child[i]];
				sumComplexKahan(divComplex(conjComplex(cinf[cr].A[j]), getComplex(cinf[cr].time[j]+cinf[cl].time[i], -x)), &ka);
				sumComplexKahan(divComplex(cinf[cr].A[j], getComplex(cinf[cr].time[j]-cinf[cl].time[i], x)), &ka);
			}
			cinf[n].A[i] = prodComplex(totalComplexKahan(&ka), cinf[cl].A[i]);
		}
		for(j=0; j<cinf[cr].sizeChild; j++) {
			int ind = j+cinf[cl].sizeChild;
			cinf[n].child[ind] = cinf[cr].child[j];
			cinf[n].time[ind] = cinf[cr].time[j]+param*tree->time[n];
			cinf[n].A[ind] = getComplex(0., 0.);
			TypeComplexKahan ka;
			initComplexKahan(&ka);
			for(i=0; i<cinf[cl].sizeChild; i++) {
				double x;
				x = ((double*)tree->info)[cinf[cl].child[i]]-((double*)tree->info)[cinf[cr].child[j]];
				sumComplexKahan(divComplex(conjComplex(cinf[cl].A[i]), getComplex(cinf[cl].time[i]+cinf[cr].time[j], -x)), &ka);
				sumComplexKahan(divComplex(cinf[cl].A[i], getComplex(cinf[cl].time[i]-cinf[cr].time[j], x)), &ka);
			}
			cinf[n].A[ind] = prodComplex(totalComplexKahan(&ka), cinf[cr].A[j]);
		}
	}
}

double getCauchyLogDensityStem(TypeCauchyInfo cinf, double *val, double start) {
	int i, *signA;
	double res, *allLogA;
	signA = (int*) malloc(cinf.sizeChild*sizeof(int));
	allLogA = (double*) malloc(cinf.sizeChild*sizeof(double));
	TypeComplexKahan ka;
	initComplexKahan(&ka);
	for(i=0; i<cinf.sizeChild; i++) {
		double costheta;
		TypeComplex logd1;
		logd1 = getComplex(cinf.time[i], val[cinf.child[i]]-start);
		costheta = cos(getImagLogComplex(cinf.A[i])-getImagLogComplex(logd1));
		signA[i] = sign(costheta);
		allLogA[i] = getRealLogComplex(cinf.A[i])-getRealLogComplex(logd1)+log(fabs(costheta));
	}
	res = log(2.) - ((double)cinf.sizeChild) * log(2 * M_PI)+ sumSignedLogTableKahan(allLogA, signA, cinf.sizeChild);
	free((void*)signA);
	free((void*)allLogA);
	return res;
}

double getCauchyLogDensityNoStem(TypeCauchyInfo cinfL, TypeCauchyInfo cinfR, double *val, double start) {
	return getCauchyLogDensityStem(cinfL, val, start)+getCauchyLogDensityStem(cinfR, val, start);
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
		densRef = getCauchyLogDensityStem(cinf[tree->root], ((double*)tree->info), start);
		if(node == tree->root) {
			for(j=0; j<nVal; j++)
				dens[j] = getCauchyLogDensityStem(cinf[tree->node[node].child], ((double*)tree->info), tabVal[j])
				+getCauchyLogDensityStem(cinf[tree->node[tree->node[node].child].sibling], ((double*)tree->info), tabVal[j])
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
				dens[j] = getCauchyLogDensityStem(cinf[tree->node[node].child], ((double*)tree->info), tabVal[j])
				+getCauchyLogDensityStem(cinf[tree->node[tree->node[node].child].sibling], ((double*)tree->info), tabVal[j])
				+getCauchyLogDensityStem(cinf[treeC->root], ((double*)tree->info), tabVal[j])
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
		densRef = getCauchyLogDensityStem(cinf[rerootedTree->root], (double*)rerootedTree->info,  start);
		freeCauchyInfo(rerootedTree->root, rerootedTree, cinf);
		fillCauchyInfo(tree->root, tree, disp, cinf);
		for(i=0; i<nVal; i++)
			dens[i] = getCauchyLogDensityNoStem(cinf[tree->node[tree->root].child], cinf[tree->node[tree->node[tree->root].child].sibling], (double*) tree->info, tabVal[i])-densRef;
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
				densRef = getCauchyLogDensityStem(cinf[rerootedTree->root], ((double*)tree->info), start);
				for(j=0; j<nVal; j++)
					dens[j] = getCauchyLogDensityStem(cinf[rerootedTree->node[rerootedTree->root].child], (double*)rerootedTree->info, start-tabVal[j]) +
						getCauchyLogDensityStem(cinf[rerootedTree->node[rerootedTree->node[rerootedTree->root].child].sibling], (double*)rerootedTree->info, start-tabVal[j])+
						getCauchyLogDensityStandard(tabVal[j], disp*tree->time[node]) - densRef;
				freeCauchyInfo(rerootedTree->root, rerootedTree, cinf);
				free((void*)cinf);
				rerootedTree->info = NULL;
				freeTree(rerootedTree);
			} else {
				cinf = (TypeCauchyInfo*) malloc(tree->sizeBuf*sizeof(TypeCauchyInfo));
				fillCauchyInfo(tree->root, tree, disp, cinf);
				densRef = getCauchyLogDensityStem(cinf[tree->root], ((double*)tree->info), start);
				for(j=0; j<nVal; j++)
					dens[j] = getCauchyLogDensityStem(cinf[tree->node[tree->root].child], (double*)tree->info, start+tabVal[j]) +
						getCauchyLogDensityStem(cinf[tree->node[tree->node[tree->root].child].sibling], (double*)tree->info, start+tabVal[j])+
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
		densRef = getCauchyLogDensityStem(cinf[tree->root], ((double*)tree->info), start);
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
			dens[j] = getCauchyLogDensityStem(cinf[tree->root], ((double*)tree->info), start) + getCauchyLogDensityStandard(tabVal[j], disp*timeSave) - densRef;
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
				densRef = getCauchyLogDensityStem(cinf[rerootedTree->root], ((double*)rerootedTree->info), ((double*)rerootedTree->info)[rerootedTree->root]);
			else
				densRef = getCauchyLogDensityStem(cinf[rerootedTree->root], ((double*)rerootedTree->info), ((double*)tree->info)[node]);				
			freeCauchyInfo(rerootedTree->root, rerootedTree, cinf);
			if(tree->node[node].child == NOSUCH) {
				int j;
				fillCauchyInfo(other, tree, disp, cinf);
				for(j=0; j<nVal; j++)
					dens[j] = getCauchyLogDensityStem(cinf[other], ((double*)tree->info), ((double*)tree->info)[node]-tabVal[j]) + getCauchyLogDensityStandard(tabVal[j], disp*tree->time[node]) - densRef;
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
					dens[j] = getCauchyLogDensityStem(cinf[rerootedTree->root], ((double*)rerootedTree->info), ((double*)rerootedTree->info)[rerootedTree->root]) + getCauchyLogDensityStandard(tabVal[j], disp*time) - densRef;
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
