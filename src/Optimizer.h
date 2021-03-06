#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include "globals.h"
#include "Tree.h"
#include "Node.h"
#include "Branch.h"
#include <set>
using namespace std;

class Optimizer
{
	public:
		Optimizer();
		virtual ~Optimizer();

		void rearrange(Tree &tree, Options &options, vector<Tree> &bestTrees);

		unsigned int optimizeSPR(Tree &tree, double cutOff, vector<Tree> &bestTrees, int maxBestTrees);
		unsigned int optimizeNNI(Tree &tree, double cutOff, vector<Tree> &bestTrees, int maxBestTrees);

		void assessTree(Tree &tree, double cutOff, vector<Tree> &bestTrees, int maxBestTrees);

		void NNI(Branch* branch, int swap, Node* root);
		void subtreePrune(Branch *fromBranch, Node *fromParent, vector<int>& insertCandidates);
		void subtreeRegraft(Branch *fromBranch, Node *fromParent, Branch *toBranch, Node *toParent, Node *root);
};

#endif /* OPTIMIZER_H_ */
