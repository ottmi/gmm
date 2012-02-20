#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include "Tree.h"
#include "Node.h"
#include "Branch.h"
using namespace std;

class Optimizer
{
	public:
		Optimizer();
		virtual ~Optimizer();

		void rearrange(Tree &tree);

		void NNI(Branch* branch, int swap);
		void SPR(Branch *fromBranch, Node *fromParent, Branch *toBranch, Node *toParent);
};

#endif /* OPTIMIZER_H_ */