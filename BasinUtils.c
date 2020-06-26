#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "BasinDefs.h"

BasinSolution *basin_classify(
	BasinSolution *lrptr,
	BasinPeriodData *pptr,
	double eps)
{
	BasinSolution *basin_newsolution();
	static int id = 0;

	if (lrptr == NULL){
		lrptr = basin_newsolution();
		lrptr->id = id;
		pptr->id = id;
		lrptr->period = pptr->period;
		lrptr->freq = 1;
		lrptr->data = pptr->data;
		lrptr->left = lrptr->right = lrptr->equal = NULL;
		id++;
	} else {
		if  (pptr->period == lrptr->period){
			if (basin_comp_period(pptr->data, lrptr->data,
				pptr->period, eps)){
				lrptr->freq++;
				pptr->id = lrptr->id;
				free(pptr->data);
			} else {
				lrptr->equal = basin_classify(lrptr->equal, pptr, eps);
			}
		} else if (pptr->period < lrptr->period){
			lrptr->left = basin_classify(lrptr->left, pptr, eps);
		} else if (pptr->period > lrptr->period){
			lrptr->right = basin_classify(lrptr->right, pptr, eps);
		}
	}
	return lrptr;
}

int basin_comp_period( BasinPoint *c1, BasinPoint *c2, int period, double eps)
{
	int i;
	int basin_comp(BasinPoint *, BasinPoint *, double);

	for (i = 0; i < period; i++){
		if (!basin_comp(c1++, c2++, eps)){
			return 0;
		}
	}
	return 1;
}

int basin_comp(BasinPoint *c1, BasinPoint *c2, double eps)
{ 
	return fabs(c1->x - c2->x) + fabs(c1->y - c2->y) < eps ? 1 : 0; 
}

void printperiod( BasinSolution *p, BasinClassifiedList *c)
{
	int i;
	if (p != NULL){
		printperiod(p->left, c);
		printperiod(p->equal, c);
		printf("id %d, period %d, freq %d\n",
			p->id, p->period, p->freq);
		for (i = 0; i < p->period; i++){
			printf("\t(%lf %lf)\n", p->data->x, p->data->y);
			p->data++;
		}
		printperiod(p->right, c);
		(c + p->id)->id = p->id;
		(c + p->id)->period = p->period;
		(c + p->id)->freq = p->freq;
	}
}

int chknum(BasinSolution *tree)
{
	static int i = 0;
	if (tree != NULL){
		chknum(tree->left);
		chknum(tree->equal);
		if (tree->id > i ){ i = tree->id;}
		chknum(tree->right);
	}
	return i;
}

void sort(BasinPoint *p, int period)
{
	int i, j, k;
	double x0, y0;

	/* BAKA sort */

	for (i = 0; i < period - 1; i++){
		for (j = i + 1; j <  period ; j++){
			if ((p+i)->x > (p+j)->x){
			x0 = (p+j)->x; y0 = (p+j)->y;
			(p+j)->x = (p+i)->x; (p+j)->y = (p+i)->y;
			(p+i)->x = x0; (p+i)->y = y0;
			}
		}
	}
}

BasinPoint *basin_newpoints(int n)
{
	BasinPoint *ret;
	if ((ret = (BasinPoint *)calloc(n,
		sizeof(BasinPoint))) == NULL){
			fprintf(stderr, 
				"cannot allocate memory for %d points\n", n);
			exit(0);
	}
	return ret;
}

BasinSolution *basin_newsolution(void)
{
	BasinSolution *ptr;

	if ((ptr = (BasinSolution *)calloc(1,
		sizeof(BasinSolution)))== NULL){
			fprintf(stderr, 
				"cannot allocate memory for newsolution\n");
			exit(0);
	}
	return ptr;
}


