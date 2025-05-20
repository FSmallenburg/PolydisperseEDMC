// Event driven MD, headers.

typedef struct sevent
{
	double time;
	struct sevent* child1;
	struct sevent* child2;
	struct sevent* parent;
	struct sparticle* p1;
	struct sparticle* p2;
	struct sevent* prevq, * nextq;
	uint8_t type;
	int queue;
	unsigned int counter2;
} event;


typedef struct sparticle
{
	double x, y, z;
	double vx, vy, vz;
	double xn, yn, zn;
	double r, vr;
	//struct sneigh* neighbors[MAXNEIGH];
	struct sneigh* neighbors;
	uint8_t nneigh;
	uint8_t nnneigh;
	double t;
	double rtarget;
	double mass;
	uint8_t edge;
	uint8_t cellx, celly, cellz;
	int boxestraveledx, boxestraveledy, boxestraveledz;
	event* firstcollision;
	unsigned int counter;
	struct sparticle* prev, * next;
	int number;
	uint8_t type;
	int binnumber;
} particle;


typedef struct sneigh
{
	struct sparticle* part;
	double dist;
} neigh;


typedef struct cmplx
{
	double re;
	double im;
} cmplx;


int main();
void printstuff();
void init();

int mygetline(char* str, FILE* f);

void randomparticles();
void randommovement();
void update(particle*);
void initcelllist();
void removefromcelllist(particle*);
void addtocelllist(particle*);

void step();
void bumpradius(particle* p1);

void makeneighborlist(particle* p1, int firsttime);
double findneighborlistupdate(particle* p1);

double findcollision(particle*, particle*, double);
double findradialcollision(particle* p1, double tmin);
void findcollisions(particle*);
void findallcollisions();
void collision(event*);

void initeventpool();
void addeventtotree(event* newevent);
void addevent(event* newevent);
event* createevent(double time, particle* p1, particle* p2, int type);
void addnexteventlist();
void removeevent(event*);

void showtree();
void shownode(event*);
void checktree();
int checknode(event*);
int overlap(particle*);
int overlaplist(particle* part, int error);
void write();
void checkcells();
void backinbox(particle* p);
void thermostat(event* ev);
void writelast();
double random_gaussian();

double find_locmin_p3();
double find_realroots_p3(double mu, int k);
double findradialcollision_tailored(particle* p1, double tmin);

int parse_input(int argc, char* argv[]);
void fccparticles();
void fccaparticles_coex();
void fccbparticles_coex();
void loadparticles();

void kickback(event* ev);

double dotprodsum(cmplx *vec1, cmplx *vec2, int l);
double minpow(int m);
double gammln(double xx);
double facs(int l, int m);
float plgndr(int l, int m, double x);
void calc_harmonics(int l, int n1, int n2, cmplx *res1, cmplx *res2);
void calc_order();
void calc_conn();
void buildsann();
