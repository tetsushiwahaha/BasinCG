#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <stdio.h>
#include <png.h>
#include "BasinDefs.h"

#define TYPE 160
#define MAXPERIOD 400
#define RED	0x0000ff
#define GREEN 0x00ff00
#define BLUE 0xff0000

BasinColList coltab[] = {
	240.0,		1.0,	1,		/* blue */
	0.0,		1.0,	1,		/* red */
	300.0,		1.0,	1,		/* magenta */
	120.0,		1.0,	1,		/* green */
	180.0,		1.0,	1,		/* cyan */
	60.0,		1.0,	1,		/* yellow */
	240.0,		0.0,	1,		/* white */
	248.0,		0.56,	1,		/* slate blue */
	350.0,		0.25,	1,		/* pink */
	276.0,		0.86,	1,		/* purple */
	83.0,		0.81,	1,		/* yellowgreen */
	203.0,		0.46,	1,		/* sky blue */
	39.0,		1.0,	1,		/* orange */
}; 

int save_png(char *, BasinHead *, BasinUCRGBList *);
void hsi2rgb(BasinDBLRGBList *rgb, double H, double S, double v);


int main(int argc, char **argv)
{	
	FILE *fopen(), *fpchunk, *fpz;
	int fd;
	char basename[BUFSIZ] = "/tmp/BasinPutXXXXXX";
	char *str;
	char fname[BUFSIZ];
    unsigned int depth, width, height;
	double bias[TYPE];
	double x, y, dx, dy;
	double in, hue, sat, hdelta, satdelta;
	double chaoshue, chaossat, divhue, divsat;
	double chaosin, divin;
	double quasihue, quasisat, quasiin;
	double arg, autob;
	int i, j, k, ix, iy;
	int autobias;
	int min[TYPE], max[TYPE], up[TYPE], low[TYPE];
	int pflag[TYPE], format, screen, test, key;
	int autoflag, biasflag, colclass, initflag[TYPE];
	int delta[TYPE], beta, sum, maxperiod;
	BasinUCRGBList *data, *dptr, *data0, *dptr0; 

	BasinColList *color, *colptr;
	BasinSolutionList *slist, *sptr;
	BasinHead *h;
	BasinClassifiedList *clist, *cptr;
	BasinDBLRGBList RGB;

	if (argc != 2){
		fprintf(stderr, "usage: %s datafile[.baz]\n", argv[0]);
		exit(1);
	}
	strcpy(fname, argv[1]);
	str = rindex(fname, '.');
	if (str == NULL){
		strcat(fname, ".baz");
	} else {
		str++;
		if (!strcmp(str, ".baz")){
			fprintf(stderr, "%s is not baz file\n", argv[1]);
			exit(-1);
		}
	}
		
	if((fpz = fopen(fname, "r"))==NULL){
		perror("fopen ");
		exit(-1);
	}
    if((fd = mkstemp(basename))== -1){
		perror("mkstemp ");
		exit(-1);
	}

	if((fpchunk = fdopen(fd, "w")) == NULL){
		perror("fdopen ");
		exit(-1);
	}

	if ((h = (BasinHead *)calloc(1, sizeof (BasinHead))) == NULL){
		perror("calloc ");
		exit(-1);
	}

	read_zfile(fpz, fpchunk);

	fclose(fpchunk);

	if((fpchunk = fopen(basename, "r"))==NULL){
		perror("fopen");
		exit(-1);
	}

	unlink(basename);

	fread(h, sizeof(BasinHead), 1, fpchunk);

	printf("header length = %d\n", h->hlength);
	printf("%d types of solutions\n", h->num);
	printf("ix = %d\n", h->ix);
	printf("iy = %d\n", h->iy);

	if (h->id != BASIN_MAGIC){ 
		fprintf(stderr, "not basin data file.\n");
		exit(-1);
	}

	if ((clist = (BasinClassifiedList *)calloc(h->num,
			sizeof (BasinClassifiedList))) == NULL){
		fprintf(stderr, "cannot allocate memory for clist\n");
		exit(-1);
	}
	printf("clist = %p\n", clist);

	fread(clist, sizeof(BasinClassifiedList), h->num, fpchunk);
	cptr = clist;

	if ((slist = (BasinSolutionList *)calloc(h->ix * h->iy,
			sizeof (BasinSolutionList))) == NULL){
		fprintf(stderr, "cannot allocate memory for sollist\n");
		exit(-1);
	}
	fread(slist, sizeof(BasinSolutionList), h->ix * h->iy, fpchunk);

	sptr = slist;

	if ((data = (BasinUCRGBList *)calloc(h->ix * h->iy, 
		sizeof (BasinUCRGBList))) == NULL){
		fprintf(stderr, "cannot allocate memory for data\n");
		exit(-1);
	}
	dptr = data;
	if ((data0 = (BasinUCRGBList *)calloc(h->ix * h->iy, 
		sizeof (BasinUCRGBList))) == NULL){
		fprintf(stderr, "cannot allocate memory\n");
		exit(-1);
	}
	dptr0 = data0;

	for (i = 0; i < h->num; i++){
	/*
		printf("id:%d period:%d freq:%d\n", 
			cptr->id, cptr->period, cptr->freq);
	*/
		cptr++;
	}
	cptr = clist;
	for(i=0; i < h->num; i++){ initflag[i] = 0;}

	sum = 0; maxperiod = 0;
	for(i = 0; i < h->iy; i++){
		for(j = 0; j < h->ix; j++){
			if (!initflag[sptr->id]){
				min[sptr->id] = max[sptr->id] = sptr->converge;
				initflag[sptr->id] = 1;
				if (sptr->period > maxperiod) maxperiod = sptr->period;
				sptr++; sum ++;
				continue;
			} else {
				if (min[sptr->id] > sptr->converge) 
					min[sptr->id] = sptr->converge;
				if (max[sptr->id] < sptr->converge)
					max[sptr->id] = sptr->converge;
				if (sptr->period > maxperiod) maxperiod = sptr->period;
				/*
				printf("sptr->id = %d\n", sptr->id);
				*/
				sptr++; sum ++;
			}
		}
	}
	if (maxperiod == 0) maxperiod = 1;
	for(i=0; i < h->num; i++){ 
		delta[i] = (max[i] - min[i]) * (max[i] - min[i]);
		if (delta[i] == 0) delta[i] = 1;
		printf("delta[%d] = %d\n", i, delta[i]);
	}
	if ((color = (BasinColList *)calloc(h->num * maxperiod, 
		sizeof (BasinColList))) == NULL){
		fprintf(stderr, "cannot allocate memory\n");
		exit(-1);
	}
	colptr = color;

	sptr = slist;

	printf("ignore the depth ? (y:1, n:0) :");
	scanf("%d", &beta);
	printf("enable automatic color assigning ? (y:1, n:0): ");
	scanf("%d", &autoflag);
	printf("using individual color for every period ? (y:1, n:0):  ");
	scanf("%d", &colclass);
	printf("add automatic bias ? (y:1, n:0):");
	scanf("%d", &autobias);
	if (autobias){
		printf("bias = ");
		scanf("%lf", &autob);
	}
	printf("add bias ? (y:1, n:0) ");
	scanf("%d", &biasflag);
	for (i = 0; i < h->num; i++) bias[i] = 0;
	if (biasflag){
		for(i=0; i < h->num; i++){ 
			int ret;
			printf("period:%d\n", (clist + i)->period); 
			printf("min[%d] = %d\n", i, min[i]); 
			printf("max[%d] = %d\n", i, max[i]); 
			printf("bias = "); 
			ret = scanf("%lf", &bias[i]); 
			if (ret != 1){ getchar();}
		}
	}

	if (!autoflag){
		for (i = 0; i < h->num; i++){
			if (colclass && cptr->period > 0){
				for (j = 0; j < abs(cptr->period); j++){
					printf("period %d: id:%d-%d, f(%d): h s i=",
						cptr->period, i, j, cptr->freq);
					scanf("%lf %lf %d", 
						&((color + i * maxperiod + j)->hue),
						&((color + i * maxperiod + j)->sat),
						&((color + i * maxperiod + j)->invert));
				} 
			} else {
				double a, b;
				int l;
					printf("period %d id:%d f(%d): h s i =",
						cptr->period, i, cptr->freq);
					scanf("%lf %lf %d", &a, &b, &l);

					for (j = 0; j < abs(cptr->period); j++){
						(color + i * maxperiod + j)->hue = a;
						(color + i * maxperiod + j)->sat = b;
						(color + i * maxperiod + j)->invert = l;
					}
			}
			cptr++;
		}
	} else {
		printf("-1 : hue & sat & in = ");
		scanf("%lf %lf %lf", &quasihue, &quasisat, &quasiin);
		printf("-2 : hue & sat & in = ");
		scanf("%lf %lf %lf", &chaoshue, &chaossat, &chaosin);
		printf("-3 : hue & sat & in = ");
		scanf("%lf %lf %lf", &divhue, &divsat, &divin);
	}

	cptr = clist;

	for(i = 0; i < h->iy; i++){
		for(j = 0; j < h->ix; j++){
			register int a, a0, a1, a2, aa0, aa1;

			in = (double)((sptr->converge - min[sptr->id])
			 * (sptr->converge - min[sptr->id] ))
				/ delta[sptr->id];
			if (autoflag) {
				switch (sptr->period){
					case -1:
						hue = chaoshue;
						sat = chaossat;
						in = chaosin;
						break;
					case -2:
						hue = divhue;
						sat = divsat;
						in = divin;
						break;
					case -3:
						hue = quasihue;
						sat = quasisat;
						in = quasiin;
						break;
					default:
						sat = coltab[(sptr->period - 1) % 13].sat;
						hue = coltab[(sptr->period - 1) % 13].hue;
						in = 1.0 - in;
						break;
				}
			} else {
				sat = (color + sptr->id * maxperiod + sptr->type)->sat;
				hue = (color + sptr->id * maxperiod + sptr->type)->hue;
				if ((color + sptr->id 
					* maxperiod + sptr->type)->invert){ in = 1.0 - in; }

			}
			if (biasflag){
				double inn;
				/*
				inn = in * in;
				in = inn;
				*/
				in *= bias[sptr->id]; if (in > 1.0) in = 1.0;
			}
			if (autobias && sptr->period > 0){ 
				in *= autob; if (in > 1.0) in = 1.0;
			}
			if (min[sptr->id] == max[sptr->id] && sptr->period > 0)
				in = 1.0;
			if (beta){ in = 1.0;}

			hsi2rgb(&RGB, hue, sat, in);

			{ 
				register int r, g, b;

				r = (int)(RGB.r * 255.0);
				g = (int)(RGB.g * 255.0);
				b = (int)(RGB.b * 255.0);

				if (r > 255) r = 255;
				if (g > 255) g = 255;
				if (b > 255) b = 255;

				dptr->r = (unsigned char)r;
				dptr->g = (unsigned char)g;
				dptr->b = (unsigned char)b;

				dptr++; 
				sptr++;
			}
		}
	}
	dptr = data;
	
	/* fit to display coordinate system */
	/* (transpose) */
	for(i = h->iy - 1; i >= 0; i--){
		for(j = 0; j < h->ix; j++){
		dptr0->r = (dptr + i * h->ix + j)->r;
		dptr0->g = (dptr + i * h->ix + j)->g;
		dptr0->b = (dptr + i * h->ix + j)->b;
		dptr0++;
		}
	}

	printf("\nsaving pngfile...");

	save_png(fname, h, data0);

	printf("done.\n");

	free(data);
	free(data0);
	free(slist); 
	free(clist); 

	return 0;
}


int save_png(char *fname, BasinHead *h, BasinUCRGBList *data)
{
	FILE *fpout;
	char buf[BUFSIZ];
	char *sptr;
	unsigned char **image;
	png_structp	png_ptr;
	png_infop	info_ptr;

	int i, j;
	int s;

	strcpy(buf, fname);
	if ((sptr = rindex(buf, '.')) != NULL){ *sptr = '\0'; }
	strcat(buf, ".png");

	if ((fpout = fopen(buf,"w"))==NULL){
		fprintf(stderr, "cannot open %s\n", buf);
		exit(0);
	}

	image = (png_bytepp)malloc(h->iy * sizeof (png_bytep)); 

	for (i = 0; i < h->iy; i++){
		image[i] = (png_bytep)malloc(h->ix * sizeof (png_byte) * 3);
		/* for RGB pixels */
	}

	for (i = 0; i < h->iy; i++){
		for (j = 0; j < h->ix; j++){
			image[i][3 * j]     = data->r;
			image[i][3 * j + 1] = data->g;
			image[i][3 * j + 2] = data->b;
			data++;
		}
	}

	png_ptr = png_create_write_struct(
		PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	info_ptr = png_create_info_struct(png_ptr);
	png_init_io(png_ptr, fpout);
	png_set_IHDR(png_ptr, info_ptr, h->ix, h->iy, 
		8 * sizeof (unsigned char), PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);
	png_write_image(png_ptr, image);
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	fclose(fpout);
	return 0;
}

void hsi2rgb(BasinDBLRGBList *rgb, double H, double S, double v)
{
	int icol;
	double h1, f;
	double a[7];

	if (H > 360.0) H = 360 - H;
	h1 = H / 60.0;
	icol = (int)h1;
	f = h1 - (double)icol;

	a[1] = v;
	a[2] = v;
	a[3] = v * (1.0 - (S * f));
	a[4] = v * (1.0 - S);
	a[5] = a[4];
	a[6] = v * (1.0 - (S * (1.0 - f)));

	if (icol > 4) icol -= 4; else icol += 2;
	rgb->r = a[icol];
	if (icol > 4) icol -= 4; else icol += 2;
	rgb->b = a[icol];
	if (icol > 4) icol -= 4; else icol += 2;
	rgb->g = a[icol];
}

