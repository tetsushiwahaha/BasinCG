#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "BasinDefs.h"

#define	CPP_COMMAND_LINE "/usr/bin/cpp -P %s"

#define MAXPERIOD 2000

int main(int argc, char **argv)
{
	char cmd[BUFSIZ], fname[BUFSIZ];
	int fd;
	char basename[BUFSIZ] = "/tmp/basin.XXXXXX";
	int mapmax, period, convflag, num;
	int i, ii, j, k, ix, iy;
	double f(), g(), x, y, x0, y0, x1, y1, dx, dy, x00, y00;
	double xmin, xmax, ymin, ymax;
	double a, b, c, eps, eps_ex;
	double delta, divergence;
	BasinPoint *data, *ptr; 
	BasinPeriodData *sol, *fixed();
	BasinSolution *root;
	BasinSolutionList *slist, *sptr;
	BasinHead *header;
	BasinClassifiedList *clist;
	FILE *fpin, *fpout;
	void sort(BasinPoint *, int);
	int chknum(BasinSolution *);
	void printperiod(BasinSolution *, BasinClassifiedList *);

	int count; 

	root = NULL;

	if (argc != 2){
		fprintf(stderr, "usage: %s datafile\n", argv[0]);
		exit(1);
	}
		
	sprintf(cmd, CPP_COMMAND_LINE, argv[1]);

	if ((fpin = popen(cmd, "r")) == NULL){
		perror("popen; ");
		exit(-1);
	}

	if ((data = (BasinPoint *)calloc(MAXPERIOD * 2,
			sizeof (BasinPoint))) == NULL){
		perror("calloc; ");
		exit(-1);
	}
	ptr = data;

	if ((sol = (BasinPeriodData *)calloc(1,
		sizeof(BasinPeriodData)))== NULL){
			perror("calloc; ");
			exit(0);
	}

	fscanf(fpin, "%lf %lf %lf", &a, &b, &c);
	fscanf(fpin, "%lf %lf", &xmin, &xmax);
	fscanf(fpin, "%lf %lf", &ymin, &ymax);
	fscanf(fpin, "%d %d", &ix, &iy);
	fscanf(fpin, "%d %lf", &mapmax, &divergence);
	fscanf(fpin, "%lf", &eps);
	fscanf(fpin, "%lf", &eps_ex);
	pclose(fpin);

	if((fd = mkstemp(basename)) == -1){
		perror("mkstemp; ");
		exit(-1);
	}

	if((fpout = fdopen(fd, "w")) == NULL){
		perror("fdopen ");
		exit(-1);
	}

	dx = (xmax - xmin) / ix;
	dy = (ymax - ymin) / iy;

	/* make a buffer to memory points */
	if ((slist = (BasinSolutionList *)calloc(ix * iy,
			sizeof (BasinSolutionList))) == NULL){
			perror("calloc; ");
		exit(-1);
	}
	sptr = slist;
	if ((header = (BasinHead *)calloc(1, sizeof (BasinHead))) == NULL){
			perror("calloc; ");
		exit(-1);
	}
	header->id = BASIN_MAGIC;
	header->ix = ix;
	header->iy = iy;

	for (i = 0; i < iy; i++){
		printf("\r(%8d rows)", i); fflush(stdout);
		for (j = 0; j < ix; j++){
			x = xmin + j * dx;
			y = ymin + i * dy;
			x00 = x; x0 = x; y00 = y; y0 = y;
			k = 0; 
			ptr = data;
			while (1) {
				double xx, yy, xx0, yy0, xx1,yy1;
		 		xx = f(x, y, a, b, c); yy = g(x, y, a, b, c);
				x = xx; y = yy;
		 		xx0 = f(x0, y0, a, b, c); yy0 = g(x0, y0, a, b, c);
				x1 = xx0; y1 = yy0;
		 		xx0 = f(x1, y1, a, b, c); yy0 = g(x1, y1, a, b, c);
				x0 = xx0; y0 = yy0;
				delta = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
				if (delta < eps){
					convflag = 0; break;
				} else if (delta > divergence){
					convflag = 2; break;
				}
				if (k > mapmax){ convflag = 1; break; }
				k++;
			}
			if (convflag == 1){
#ifdef MSG
				printf("chaotic solution found\n");
#endif
				sol->data = NULL;
				sol->period = -1;
				root = basin_classify(root, sol, eps_ex);
				sptr->period = sol->period;
				sptr->id = sol->id;
				sptr->converge = 0;
				sptr++;
			} else if (convflag == 2){
#ifdef MSG
				printf("divergence!\n");
#endif
				sol->data = NULL;
				sol->period = -2;
				root = basin_classify(root, sol, eps_ex);
				sptr->period = sol->period;
				sptr->id = sol->id;
				sptr->converge = 0;
				sptr++;
			} else {
				double xx, yy;
				x0 = x; y0 = y;
				period = 0;
				do {
		 			xx = f(x, y, a, b, c);
					yy = g(x, y, a, b, c);
					x = xx; y = yy;
					ptr->x = x; ptr->y = y;
					ptr++;
					period++;
					delta = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
				} while (delta > eps_ex);
				ptr = data;
				sol->period = period;
				sol->data = basin_newpoints(period);
				for (ii = 0; ii < period; ii++){
					(sol->data + ii)->x = (ptr + ii)->x;
					(sol->data + ii)->y = (ptr + ii)->y;
				}
				sort(sol->data, period);
				root = basin_classify(root, sol, eps_ex);
				x = x00; x0 = x00; y = y00; y0 = y00;
				while (1) {
					double xx, yy, xx0, yy0, xx1, yy1;
					int ii;
					for (ii = 0; ii < period; ii++){
			 			xx = f(x, y, a, b, c);
						yy = g(x, y, a, b, c);
						x = xx; y = yy;
			 			xx0 = f(x0, y0, a, b, c);
						yy0 = g(x0, y0, a, b, c);
						x1 = xx0; y1 = yy0;
			 			xx0 = f(x1, y1, a, b, c);
						yy0 = g(x1, y1, a, b, c);
						x0 = xx0; y0 = yy0;
						delta = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
					}
					if (delta < eps_ex){
						break;
					}
				}
				for (ii = 0; ii < period; ii++){
					if (fabs((sol->data + ii)->x - x) + 
						fabs((sol->data + ii)->y - y) < eps_ex * 2.0){
						sptr->type = ii;
					}
				}
				sptr->period = sol->period;
				sptr->id = sol->id;
				sptr->converge = k;
				sptr++;
			}
		}
	}	
	header->num = chknum(root) + 1;
	printf("found %d periodic solutions\n", header->num);
	if ((clist = (BasinClassifiedList *)calloc(header->num,
			sizeof (BasinClassifiedList))) == NULL){
		fprintf(stderr, 
			"cannot allocate memory for classification list\n");
		exit(-1);
	}
	header->hlength = sizeof (BasinClassifiedList) * (header->num);

	fwrite(header, sizeof(BasinHead), 1, fpout);
	printperiod(root, clist);
	fwrite(clist, sizeof(BasinClassifiedList), header->num, fpout);
	fwrite(slist, sizeof(BasinSolutionList), ix * iy, fpout);
	fclose(fpout);

	if((fpin = fopen(basename, "r"))==NULL){
		perror("fopen ");
		exit(-1);
	}
	unlink(basename);

	sprintf(fname, "%s.baz", argv[1]);
	if((fpout = fopen(fname, "w")) == NULL){
		fprintf(stderr, "cannot open %s\n",argv[1]);
		exit(-1);
	}
	write_zfile(fpin, fpout);
	fclose(fpin);
	fclose(fpout);
}

