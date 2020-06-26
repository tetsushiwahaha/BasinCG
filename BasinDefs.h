typedef struct {
	int id;
	int hlength;
	int num;
	int ix;
	int iy;
} BasinHead;

typedef struct {
	double x;
	double y;
} BasinPoint;

typedef struct {
	int reason;
	int period;
	int id;
	int converge;
	BasinPoint *data;
} BasinPeriodData;

typedef struct _basinSolution {
	int id;
	int freq;
	int period;
	int color;
	int maxconv;
	BasinPoint *data;
	struct _basinSolution *left;
	struct _basinSolution *right;
	struct _basinSolution *equal;
} BasinSolution;

typedef struct {
	int id;
	int type;
	int period;
	int converge;
} BasinSolutionList;

typedef struct {
	int id;
	int period;
	int freq;
} BasinClassifiedList;

typedef struct {
	double r;
	double g;
	double b;
} BasinDBLRGBList;

typedef struct {
	unsigned char r;
	unsigned char g;
	unsigned char b;
} BasinUCRGBList;

typedef struct {
	double h;
	double s;
	double i;
} BasinHSIList;

typedef struct {
	double hue;
	double sat;
	int invert;
} BasinColList;

#define BASIN_MAGIC 0xaabbccdd

void read_zfile(FILE *, FILE *);
void write_zfile(FILE *, FILE *);

BasinSolution *basin_classify( BasinSolution *, BasinPeriodData *, double eps);
int basin_comp_period(BasinPoint *, BasinPoint *, int, double);

BasinPoint *basin_newpoints(int n);

