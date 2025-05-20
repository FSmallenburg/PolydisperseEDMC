#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>


//Include random number generator 
#include "xoshiro256plus.c"

//Include algorithm for exponential random numbers 
#include "ran_exp.c"


#include "EDMC_poly.h"


//Size of array allocated for the event tree (currently overkill)
#define MAXEVENTS (2*N+20)

//Maximum number of cells in each direction
#define CEL 75

//Pi (if not already defined)
#ifndef M_PI
#define M_PI 3.1415926535897932
#endif
#define INVFPI 0.07957747154594766788444188168625718101 // = 1 / (4 * pi)

//Number of particles:
#define N 16384
//Actual number of particles in use
int Np = N;

double maxtime = 100;
double density = 0.9;               


//We need a maximum particle size so that we can make a proper cell list. This should be large enough compared to the size distribution of the system.
//It is probably also good to have a minimum size, so particles cannot hit size zero.
//Both these values are used in predicting the interactions of the radii with the field.
double maxparticleradius = 0.75;
double minparticleradius = 0.05;

int makesnapshots = 1;         //Whether to make snapshots during the run (yes = 1, no = 0)
double writeinterval = 1;      //Time between output to screen / data file
double snapshotinterval = 10;  //Time between snapshots (should be a multiple of writeinterval)


//Variables related to the event queueing system. These can affect efficiency.
//The system schedules only events in the current block of time with length "eventlisttime" into a sorted binary search tree. 
//The rest are scheduled in unordered linked lists associated with the "numeventlists" next blocks.
//"numeventlists" is roughly equal to maxscheduletime / eventlisttime
//Any events occurring even later are put into an overflow list
//After every time block with length "eventlisttime", the set of events in the next linear list is moved into the binary search try.
//All events in the overflow list are also rescheduled.

//After every "writeinterval", the code will output two listsizes to screen. 
//The first is the average number of events in the first that gets moved into the event tree after each block.
//The second is the average length of the overflow list.
//Ideally, we set maxscheduletime large enough that the average overflow list size is negligible (i.e. <10 events)
//Also, there is some optimum value for the number of events per block (scales approximately linearly with "eventlisttime").
//I seem to get good results with an eventlisttime chosen such that there are a few hundred events per block, and dependence is pretty weak (similar performance in the range of e.g. 5 to 500 events per block...)
#define MAXNUMEVENTLISTS (10*N)
double maxscheduletime = 5;
int numeventlists;
double eventlisttime = 5.0f / (float)N;


//Thermostat
int usethermostat = 0;              //Whether to use a thermostat
double thermostatinterval = 0.01;   //Time between applications of thermostat, which gets rid of excess heat generated while growing



//Size of neighbor lists
double shellsize = 1.2;
//Maximum number of neighbors in neighbor list. Simulation will crash if the number of particles in a neighbor list exceeds this number
#define MAXNEIGH 100


//Internal variables
double time = 0;
double reftime = 0;
int currentlist = 0;
const double never = 9999999999999999999999.9;
int listcounter1 = 0, listcounter2 = 0, mergecounter = 0;
int nbig;

event* eventlists[MAXNUMEVENTLISTS + 1];    //Last one is overflow list

particle particles[N];
particle particleseq[N];
particle* celllist[CEL][CEL][CEL];
event eventlist[MAXEVENTS];
event* root;
event* eventpool[MAXEVENTS];
int nempty = 0;
double xsize, ysize, zsize;     //Box size
double hx, hy, hz;              //Half box size
double cxsize, cysize, czsize;  //Cell size
int    cx, cy, cz;              //Number of cells
double dvtot = 0;               //Momentum transfer (for calculating pressure)
unsigned int colcounter = 0;    //Collision counter (will probably overflow in a long run...)
int stop = 0;


//Kickback -- to pin coex. interface in place
int usekickback = 0;            //Whether to use kickback
double kickinterval = 0.05;
double binsize = 0.0400;
double Xoffset = 0.0f;
char binfilename[200];
FILE* binfile = NULL;


//RNG
unsigned long seed = 1;


//Radial collision with external field (chem. pot. diff. which is here a 3rd order polynomial)
// 
double a0 = 0.0;
double a1 = 0.0;
double a2 = 0.0;
double a3 = 0.0;
double rmin = 0.0;
double umin = 0.0;
double depressed_p = 0.0f;
double* roots = NULL;


//Stress tensor
double Sij[3] = {0, 0, 0};


//Choice of particle initialization
char inputfilename[255] = "\0";         //Name of input configuration file
int initchoice = 0;                     //Choice of initialization, default is random particles ('0'), see '-h' for more
double FCClatticeparam = sqrt(2.0f);    //Default value for monodisperse system with unity particle size
double fdensity = 0.0f;            //Fluid trial density
double Xdensity = 0.0f;            //Crystal trial density
double propVolX = 0.5f;                 //Proportion of box volume to fill with crystal, defaut is 1/2


//Verbose
int verbose = 0;
char vfilename[200];
FILE* vfile = NULL;


//BOPs -- ten Wolde's q6 & Lechner and Dellago's barq6
cmplx *bop6 = NULL;
cmplx *qbar = NULL;


int main(int argc, char* argv[])
{
    parse_input(argc, argv);
    init();
    if (stop) return 1;				//Stop is set to 1 whenever the simulation should stop
    printf("[EDMC-poly] Starting\n");
    while (!stop)
    {
        step();
    }
    printstuff();
    return 0;
}


/**************************************************
**                  PARSE_INPUT
**************************************************/
int parse_input(int argc, char* argv[]){
    struct option longopt[] = {
        {"density",required_argument,NULL,'d'},
        {"maxtime",required_argument,NULL,'m'},
        {"randseed",required_argument,NULL,'s'},
        {"a0",required_argument,NULL,'0'},
        {"a1",required_argument,NULL,'1'},
        {"a2",required_argument,NULL,'2'},
        {"a3",required_argument,NULL,'3'},
        {"init",required_argument,NULL,'i'},
        {"FCClatticeparam",required_argument,NULL,'a'},
        {"fdensity",required_argument,NULL,'f'},
        {"Xdensity",required_argument,NULL,'X'},
        {"propVolX",required_argument,NULL,'p'},
        {"input",required_argument,NULL,'I'},
        {"kickback",no_argument,NULL,'k'},
        {"verbose",no_argument,NULL,'v'},
        {"help",required_argument,NULL,'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "d:m:s:0:1:2:3:i:a:f:X:p:I:kvh", longopt, NULL)) != -1){
        switch (c){
            case 'd':
                if (sscanf(optarg, "%lf", &density) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse density value\n");
                    exit(3);
                }
                break;
            case 'm':
                if (sscanf(optarg, "%lf", &maxtime) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse maxtime value\n");
                    exit(3);
                }
                break;
            case 's':
                if (sscanf(optarg, "%lu", &seed) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse seed value\n");
                    exit(3);
                }
                break;
            case '0':
                if (sscanf(optarg, "%lf", &a0) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse a0 value\n");
                    exit(3);
                }
                break;
            case '1':
                if (sscanf(optarg, "%lf", &a1) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse a1 value\n");
                    exit(3);
                }
                break;
            case '2':
                if (sscanf(optarg, "%lf", &a2) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse a2 value\n");
                    exit(3);
                }
                break;
            case '3':
                if (sscanf(optarg, "%lf", &a3) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse a3 value\n");
                    exit(3);
                }
                break;
            case 'i':
                if (sscanf(optarg, "%d", &initchoice) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse initchoice value\n");
                    exit(3);
                }
                break;
            case 'a':
                if (sscanf(optarg, "%lf", &FCClatticeparam) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse FCClatticeparam value\n");
                    exit(3);
                }
                break;
            case 'f':
                if (sscanf(optarg, "%lf", &fdensity) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse fdensity value\n");
                    exit(3);
                }
                break;
            case 'X':
                if (sscanf(optarg, "%lf", &Xdensity) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse Xdensity value\n");
                    exit(3);
                }
                break;
            case 'p':
                if (sscanf(optarg, "%lf", &propVolX) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse propVolX value\n");
                    exit(3);
                }
                break;
            case 'I':
                if (sscanf(optarg, "%s", inputfilename) != 1){
                    printf("[EDMC-poly] ERROR: Could not parse inputfilename value\n");
                    exit(3);
                }
                break;
            case 'k':
                usekickback = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'h':
                printf("[EDMC-poly] OPTIONS:\n");
                printf("\t-d [double] target density of the whole simulation box\n");
                printf("\t-m [int] total simulation time\n");
                printf("\t-r [int] seed for random number generator\n");
                printf("\t-0 [double] a0 external field parameter\n");
                printf("\t-1 [double] a1 external field parameter\n");
                printf("\t-2 [double] a2 external field parameter\n");
                printf("\t-3 [double] a3 external field parameter\n");
                printf("\t-i [0/1/2/3/4] choice of initial configuration from: 0-random (default), 1-FCC lattice in cubic box, 2-FCCa lattice in elongated (3:1) box, 3-FCCb lattice in elongated (3:1) box, 4-load from file\n");
                printf("\t-f [double] fluid density for direct coexistence simulations\n");
                printf("\t-X [double] crystal density for direct coexistence simulations\n");
                printf("\t-p [double] proportion of the box volume to fill with crystal (default: 0.5)\n");
                printf("\t-I [string] path to inial configuration file\n");
                printf("\t-k keep coexistence plane in place (messes up travelled boxes count)\n");
                printf("\t-v save verbose to a file\n");
                printf("\t-h help\n");
                exit(0);
        }
    }

    if (initchoice == 4 && inputfilename[0] == '\0'){
        printf("[EDMC-poly] WARNING: No supplied initial configuration file name, resorting to default ('input.osph')");
        strcpy(inputfilename, "input.osph");
    }

    return 0;
}


/**************************************************
**                 PRINTSTUFF
** Some data at the end of the simulation
**************************************************/
void printstuff()
{
    int i;
    particle* p;
    double v2tot = 0;
    double vfilled = 0;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
        vfilled += p->r * p->r * p->r * 8;
    }
    vfilled *= M_PI / 6.0;
    printf("[EDMC-poly] Average kinetic energy: %lf\n", 0.5 * v2tot / N);
    double volume = xsize * ysize * zsize;
    double dens = N / volume;
    double press = -dvtot / (3.0 * volume * time);
    double pressid = dens;
    double presstot = press + pressid;
    printf("[EDMC-poly] Total time simulated  : %lf\n", time);
    printf("[EDMC-poly] Density               : %lf\n", (double) N / volume);
    printf("[EDMC-poly] Packing fraction      : %lf\n", vfilled / volume);
    printf("[EDMC-poly] Measured pressure     : %lf + %lf = %lf\n", press, pressid, presstot);

    if (verbose && vfile!=NULL) fclose(vfile);

}


/**************************************************
**                    INIT
**************************************************/
void init()
{
    random_exponential_init();

    int i;
    printf("[EDMC-poly] Seed: %u\n", (int)seed);
    printf("[EDMC-poly] fdensity: %lf\n[EDMC-poly] Xdensity: %lf\n", fdensity, Xdensity);
    init_genrand(seed);
    initeventpool();

    //Field interaction initialization
    depressed_p = (3.0f*a3*a1 - a2*a2) / (3.0f*a3*a3);
    rmin = find_locmin_p3();
    umin = a0 + a1*rmin + a2*rmin*rmin + a3*rmin*rmin*rmin;
    printf("[EDMC-poly] FIELD: rmin: %.6lf umin: %.6lf\n", rmin, umin);
    roots = malloc(3*sizeof(*roots));

    //BOP initialization
    bop6 = malloc(N * (2 * 6 + 1) * sizeof(*bop6));
    qbar = malloc(N * (2 * 6 + 1) * sizeof(*qbar));

    for (i = 0; i < N; i++)
    {
        particle* p = particles + i;
        p->number = i;
        p->boxestraveledx = 0;
        p->boxestraveledy = 0;
        p->boxestraveledz = 0;
        p->nneigh = 0;
        p->nnneigh = 0;
        p->counter = 0;
        p->t = 0;

        p->neighbors = malloc(MAXNEIGH * sizeof(*(p->neighbors)));
    }

    switch (initchoice){
        case 0:
            randomparticles();
            break;
        case 1:
            fccparticles();
            break;
        case 2:
            fccaparticles_coex();
            break;
        case 3:
            fccbparticles_coex();
            break;
        case 4:
            loadparticles();
            break;
        default:
            randomparticles();
            break;
    }
    randommovement();
    hx = 0.5 * xsize; hy = 0.5 * ysize; hz = 0.5 * zsize;	//Values used for periodic boundary conditions



    for (i = 0; i < N; i++)
    {
        particle* p = particles + i;
        p->xn = p->x;
        p->yn = p->y;
        p->zn = p->z;
        p->binnumber = (int) (p->z / binsize);
    }

    //printf("[DEBUG] Before: %d events\n", MAXEVENTS - nempty);
    for (i = 0; i < N; i++)
    {
        makeneighborlist(particles + i, 1);
    }
    printf("[EDMC-poly] Done adding collisions: %d events\n", MAXEVENTS - nempty);


    if (usethermostat)
    {
        thermostat(NULL);
        printf("[EDMC-poly] Started thermostat: %d events\n", MAXEVENTS - nempty);
    }
    
    //Verbose file initialization
    if (verbose){
        if (usekickback) sprintf(vfilename, "verbose.n%d.v%.4lf.b%.4lf.osph", N, xsize*ysize*zsize, binsize);
        else sprintf(vfilename, "verbose.n%d.v%.4lf.osph", N, xsize * ysize * zsize);
        vfile = fopen(vfilename, "w");
        fclose(vfile);
        vfile = fopen(vfilename, "a");
    }

    //Make a copy of initial system state
    for (i=0; i<N; i++) particleseq[i] = particles[i];

    if (usekickback) {
        kickback(NULL);
        printf("[EDMC-poly] Initialized set of kickback events: %d events\n", MAXEVENTS - nempty);
        sprintf(binfilename, "bins.n%d.v%.4lf.b%.4f.dat", N, xsize*ysize*zsize, binsize);
        binfile = fopen(binfilename, "w");
        fclose(binfile);
        binfile = fopen(binfilename, "a");
    }

}


/******************************************************
**               MYGETLINE
** Reads a single line, skipping over lines
** commented out with #
******************************************************/
int mygetline(char* str, FILE* f)
{
    int comment = 1;
    while (comment)
    {
        if (!fgets(str, 255, f)) return -1;
        if (str[0] != '#') comment = 0;
    }
    return 0;
}


/**************************************************
**                    RANDOMPARTICLES
** Positions particles randomly in the box
** Particles start small
**************************************************/
void randomparticles()
{
    double vol = N / density;
    printf("[EDMC-poly] Volume: %lf\n", vol);

    xsize = cbrt(vol);
    ysize = xsize;
    zsize = ysize;
    initcelllist();
    int i;
    particle* p;
    for (i = 0; i < N; i++) //First put particles at zero
    {
        particles[i].x = 0; particles[i].y = 0; particles[i].z = 0;
    }
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->rtarget = 0.5;
        p->r = 0.5 * p->rtarget;    //Start particles off small, so it's easy to make a random configuration
        p->mass = 1;
        p->type = 0;
        p->number = i;
        do
        {
            p->x = genrand_real2() * xsize;			//Random location and speed
            p->y = genrand_real2() * ysize;
            p->z = genrand_real2() * zsize;
            p->cellx = p->x / cxsize;				//Find particle's cell
            p->celly = p->y / cysize;
            p->cellz = p->z / czsize;
        } while (overlaplist(p, 0));
        double sqm = 1.0 / sqrt(p->mass);
        p->xn = p->x; p->yn = p->y; p->zn = p->z;
        p->vx = (genrand_real2() - 0.5) / sqm;
        p->vy = (genrand_real2() - 0.5) / sqm;
        p->vz = (genrand_real2() - 0.5) / sqm;
        p->t = 0;   //r and v known at t=0
        p->next = celllist[p->cellx][p->celly][p->cellz];   //Add particle to celllist
        if (p->next) p->next->prev = p; //Link up list
        celllist[p->cellx][p->celly][p->cellz] = p;
        p->prev = NULL;
    }
}


/**************************************************
**                RANDOMMOVEMENT
**************************************************/
void randommovement()
{
    particle* p;
    double v2tot = 0, vxtot = 0, vytot = 0, vztot = 0;
    double mtot = 0;
    int i;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        double imsq = 1.0 / sqrt(p->mass);

        p->vx = imsq * random_gaussian();
        p->vy = imsq * random_gaussian();
        p->vz = imsq * random_gaussian();
        vxtot += p->mass * p->vx;			//Keep track of total v
        vytot += p->mass * p->vy;
        vztot += p->mass * p->vz;
        mtot += p->mass;

        p->vr = random_gaussian();
    }


    vxtot /= mtot; vytot /= mtot; vztot /= mtot;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx -= vxtot;					    //Make sure v_cm = 0
        p->vy -= vytot;
        p->vz -= vztot;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz) + p->vr*p->vr;
    }
    double fac = sqrt(4.0 / (v2tot / N));
    v2tot = 0;
    vxtot = vytot = vztot = 0;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx *= fac;					    //Fix energy
        p->vy *= fac;
        p->vz *= fac;
        p->vr *= fac;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
        vxtot += p->mass * p->vx;			//Keep track of total v
        vytot += p->mass * p->vy;
        vztot += p->mass * p->vz;
    }
    printf("[EDMC-poly] average v2: %lf (%lf, %lf, %lf)\n", v2tot / N, vxtot / N, vytot / N, vztot / N);
}


/**************************************************
**                UPDATE
**************************************************/
void update(particle* p1)
{
    double dt = time - p1->t;
    p1->t = time;
    p1->x += dt * p1->vx;
    p1->y += dt * p1->vy;
    p1->z += dt * p1->vz;
    p1->r += dt * p1->vr;

}


/**************************************************
**                 INITCELLLIST
**************************************************/
void initcelllist()
{
    int i, j, k;
    double mincellsize = shellsize * 2*maxparticleradius;
    cx = (int)(xsize - 0.0001) / mincellsize;	//Set number of cells
    cy = (int)(ysize - 0.0001) / mincellsize;
    cz = (int)(zsize - 0.0001) / mincellsize;
    printf("[EDMC-poly] Cells: %d, %d, %d\n", cx, cy, cz);
    if (cx >= CEL || cy >= CEL || cz >= CEL)
    {
        printf("[EDMC-poly] ERROR: Too many cells!\n");
        stop = 1; return;
    }
    cxsize = xsize / cx;						//Set cell size
    cysize = ysize / cy;
    czsize = zsize / cz;
    for (i = 0; i < CEL; i++)					//Clear celllist
        for (j = 0; j < CEL; j++)
            for (k = 0; k < CEL; k++)
            {
                celllist[i][j][k] = NULL;
            }

}


/**************************************************
**               REMOVEFROMCELLLIST
**************************************************/
void removefromcelllist(particle* p1)
{
    if (p1->prev) p1->prev->next = p1->next;    //Remove particle from celllist
    else          celllist[p1->cellx][p1->celly][p1->cellz] = p1->next;
    if (p1->next) p1->next->prev = p1->prev;
}


/**************************************************
**                    ADDTOCELLLIST
**************************************************/
void addtocelllist(particle* p)
{
    p->cellx = p->x / cxsize;				//Find particle's cell
    p->celly = p->y / cysize;
    p->cellz = p->z / czsize;
    p->next = celllist[p->cellx][p->celly][p->cellz];   //Add particle to celllist
    if (p->next) p->next->prev = p;			//Link up list
    celllist[p->cellx][p->celly][p->cellz] = p;
    p->prev = NULL;
    p->edge = (p->cellx == 0 || p->celly == 0 || p->cellz == 0 || p->cellx == cx - 1 || p->celly == cy - 1 || p->cellz == cz - 1);
}


/**************************************************
**                     STEP
**************************************************/
void step()
{
    event* ev;
    ev = root->child2;
    if (ev == NULL)
    {
        addnexteventlist();
        ev = root->child2;
    }

    while (ev->child1) ev = ev->child1;     //Find first event


    if (ev->time > maxtime)
    {
        time = maxtime;
        write();
        writelast();
        printf("[EDMC-poly] Time is up!\n");
        stop = 1;
    }

    if (ev->type == 100)
    {
        time = ev->time;
        removeevent(ev);
        write();
    }
    else if (ev->type == 200)
    {
        thermostat(ev);
    }
    else if (ev->type == 150)
    {
        time = ev->time;
        removeevent(ev);
        kickback(ev);
    }
    else if (ev->type == 8)
    {
        time = ev->time;
        removeevent(ev);
        makeneighborlist(ev->p1, 0);
    }
    else if (ev->type == 9)
    {
        time = ev->time;
        removeevent(ev);
        bumpradius(ev->p1);
    }
    else
    {
        collision(ev);
    }
}


/**************************************************
**                BUMPRADIUS
**************************************************/
void bumpradius(particle* p1)
{
    update(p1);
    p1->vr *= -1;
    p1->counter++;
    findcollisions(p1);
}


/**************************************************
**                MAKENEIGHBORLIST
**************************************************/
void makeneighborlist(particle* p1, int firsttime)
{
    int cdx, cdy, cdz, cellx, celly, cellz;
    particle* p2;
    double dx, dy, dz, r2, rm;

    update(p1);

    if (p1->x >= xsize) { p1->x -= xsize; p1->boxestraveledx++; }
    else if (p1->x < 0) { p1->x += xsize; p1->boxestraveledx--; }
    if (p1->y >= ysize) { p1->y -= ysize; p1->boxestraveledy++; }
    else if (p1->y < 0) { p1->y += ysize; p1->boxestraveledy--; }
    if (p1->z >= zsize) { p1->z -= zsize; p1->boxestraveledz++; }
    else if (p1->z < 0) { p1->z += zsize; p1->boxestraveledz--; }
    p1->xn = p1->x;
    p1->yn = p1->y;
    p1->zn = p1->z;

    removefromcelllist(p1);
    addtocelllist(p1);

    int i, j;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i].part;
        for (j = 0; j < p2->nneigh; j++)
        {
            if (p2->neighbors[j].part == p1)
            {
                p2->nneigh--;
                p2->neighbors[j].part = p2->neighbors[p2->nneigh].part;
                break;
            }
        }
    }

    cellx = p1->cellx + cx;
    celly = p1->celly + cy;
    cellz = p1->cellz + cz;

    p1->nneigh = 0;

    for (cdx = cellx - 1; cdx < cellx + 2; cdx++)
        for (cdy = celly - 1; cdy < celly + 2; cdy++)
            for (cdz = cellz - 1; cdz < cellz + 2; cdz++)
            {
                p2 = celllist[cdx % cx][cdy % cy][cdz % cz];
                while (p2)
                {
                    if (p2 != p1)
                    {
                        update(p2);
                        dx = p1->xn - p2->xn;
                        dy = p1->yn - p2->yn;
                        dz = p1->zn - p2->zn;
                        if (p1->edge)
                        {
                            if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
                            if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
                            if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
                        }
                        r2 = dx * dx + dy * dy + dz * dz;
                        rm = 2*maxparticleradius * shellsize;
                        if (r2 < rm * rm)
                        {
                            if (p1->nneigh == MAXNEIGH || p2->nneigh == MAXNEIGH)
                            {
                                printf ("Too many neighbors\n");
                                exit(3);
                            }
                            p1->neighbors[p1->nneigh++].part = p2;
                            p2->neighbors[p2->nneigh++].part = p1;
                        }
                    }
                    p2 = p2->next;
                }
            }

    findcollisions(p1);
}


/**************************************************
**                FINDNEIGHBORLISTUPDATE
** Assumes p1 is up to date
** Note that the particle is always in the same
** box as its neighborlist position (p->xn)
**************************************************/
double findneighborlistupdate(particle* p1)
{
    double dx = p1->x - p1->xn;
    double dy = p1->y - p1->yn;
    double dz = p1->z - p1->zn;

    double dvx = p1->vx, dvy = p1->vy, dvz = p1->vz;

    double b = dx * dvx + dy * dvy + dz * dvz;                  //dr.dv

    double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
    double dr2 = dx * dx + dy * dy + dz * dz;
    double md = (shellsize - 1) * maxparticleradius;
    double disc = b * b - dv2 * (dr2 - md * md);
    double t = (-b + sqrt(disc)) / dv2;
    return t;
}


/**************************************************
**                FINDCOLLISION
** Detect the next collision for two particles
** Note that p1 is always up to date in
** findcollision
**************************************************/
double findcollision(particle* p1, particle* p2, double tmin)
{
    double dt2 = time - p2->t;
    double dx = p1->x - p2->x - dt2 * p2->vx;    //relative distance at current time
    double dy = p1->y - p2->y - dt2 * p2->vy;
    double dz = p1->z - p2->z - dt2 * p2->vz;
    if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
    if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
    if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
    double dvx = p1->vx - p2->vx;                               //relative velocity
    double dvy = p1->vy - p2->vy;
    double dvz = p1->vz - p2->vz;
    double dvr = p1->vr + p2->vr;
    double md = p1->r + p2->r + dt2 * p2->vr;


    double b = dx * dvx + dy * dvy + dz * dvz - dvr * md;       //dr.dv
    double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
    double twoa = dv2 - dvr * dvr;
   
    if (b > 0 && twoa > 0) return never;

    double dr2 = dx * dx + dy * dy + dz * dz;
    double disc = b * b - twoa * (dr2 - md * md);
    if (disc < 0) return never;
    double t = (-b - sqrt(disc)) / twoa;
    return t;
}


// /**************************************************
// **                FINDRADIALCOLLISION
// ** Field interaction
// ** Based on a simple quadratic external field V(R_i)
// **************************************************/
// double findradialcollision(particle* p1, double tmin)
// {
//     double rmin = 0.5;      // well minimum location
//     double kappa = 500;     // spring constant (very high = doesn't wiggle much)

//     double rcur = p1->r;
//     double vrcur = p1->vr;

//     if(vrcur > 0)   //Particle is currently growing
//     {
//         double ucurrent = 0;
//         if(rcur > rmin) ucurrent = kappa*(rcur-rmin)*(rcur-rmin);
//         double unew = ucurrent + random_exponential();      //Energy at collision
//         double rnew = rmin + sqrt(unew / kappa);
//         if (rnew > maxparticleradius) rnew = maxparticleradius;
//         double t = (rnew - rcur) / vrcur;
//         return t;
//     }
//     else            //Particle is currently shrinking
//     {
//         double ucurrent = 0;
//         if(rcur < rmin) ucurrent = kappa*(rcur-rmin)*(rcur-rmin);
//         double unew = ucurrent + random_exponential();      //Energy at collision
//         double rnew = rmin - sqrt(unew / kappa) ;
//         if (rnew < minparticleradius) rnew = minparticleradius;
//         double t = (rnew - rcur) / vrcur;
//         return t;        
//     }
// }


/**************************************************
**                FINDCOLLISIONS
** Find all collisions for particle p1.
** The particle 'not' isn't checked.
**************************************************/
void findcollisions(particle* p1)    //All collisions of particle p1
{
    int i;
    double t;
    particle* partner = p1;
    particle* p2;

    double tmin = findneighborlistupdate(p1);
    int type = 8;

    t = findradialcollision_tailored(p1, tmin);
    if (t < tmin)
    {
        tmin = t;
        type = 9;
    }
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i].part;
        t = findcollision(p1, p2, tmin);
        if (t < tmin)
        {
            tmin = t;
            partner = p2;
            type = 0;
        }
    }
    if (tmin < 0) printf("[EDMC-poly] WARNING: type: %d tmin: %lf time: %lf p1: %d p2:%d\n", type, tmin , time, p1->number, partner->number);
    if (partner->number == p1->number && type == 0) printf("[EDMC-poly] WARNING p:%d colliding with itself, tmin: %lf type: %d\n", p1->number, tmin, type);
    event* ev = createevent(tmin + time, p1, partner, type);
    p1->firstcollision = ev;
    ev->counter2 = partner->counter;
}



/**************************************************
**                FINDALLCOLLISION
** All collisions of all particle pairs
**************************************************/
void findallcollisions()       //All collisions of all particle pairs
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        particle* p1 = particles + i;
        particle* partner = p1;
        double tmin = findneighborlistupdate(p1);
        int type = 8;
        for (j = 0; j < p1->nneigh; j++)
        {
            particle* p2 = p1->neighbors[j].part;
            if (p2 > p1)
            {
                double t = findcollision(p1, p2, tmin);
                if (t < tmin)
                {
                    tmin = t;
                    partner = p2;
                    type = 0;
                }
            }
        }
        if (partner)
        {
            event* ev = createevent(tmin, p1, partner, type);
            p1->firstcollision = ev;
            ev->counter2 = partner->counter;
        }
    }
}


/**************************************************
**                  COLLISION
** Process a single collision event
**************************************************/
void collision(event* ev)
{
    time = ev->time;
    particle* p1 = ev->p1;
    particle* p2 = ev->p2;
    update(p1);
    removeevent(ev);
    if (ev->counter2 != p2->counter)
    {
        findcollisions(p1);
        return;
    }

    //     if (ev->counter1 != p1->counter)
    //     {
    //         printf("Huh?\n");
    //     }
    
    update(p2);
    p1->counter++;
    p2->counter++;

    double m1 = p1->mass, r1 = p1->r;
    double m2 = p2->mass, r2 = p2->r;

    double r = r1 + r2;
    double rinv = 1.0 / r;
    double dx = (p1->x - p2->x);			//Normalized distance vector
    double dy = (p1->y - p2->y);
    double dz = (p1->z - p2->z);
    if (p1->edge)
    {
        if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
        if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
        if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
    }
    dx *= rinv;  dy *= rinv;  dz *= rinv;

    double dvx = p1->vx - p2->vx;                               //relative velocity
    double dvy = p1->vy - p2->vy;
    double dvz = p1->vz - p2->vz;
    double dvr = p1->vr + p2->vr;

    const double mR = 1;
    double b = dx * dvx + dy * dvy + dz * dvz - dvr;            //dr.dv
    double massfactor = 2.0f / (1.0f/m1 + 1.0f/m2 + 2.0f/mR);
    b *= massfactor;
    double dv1 = b/m1;
    double dv2 = b/m2;
    dvtot += r*b;
    
    //Update stress tensor
    double fac = m1*dv1*r;
    Sij[0] += fac*dx*dx;
    Sij[1] += fac*dy*dy;
    Sij[2] += fac*dz*dz;

    if (b> 0)
    {
        printf("[EDMC-poly] time: %lf B: %lf  (%d, %d)\n", time, b, p1->number, p2->number);
        exit(3);
    } 

    p1->vx -= dv1 * dx;         //Change velocities after collision
    p1->vy -= dv1 * dy;         //delta v = (-) dx2.dv2
    p1->vz -= dv1 * dz;
    p1->vr += b/mR;
    p2->vx += dv2 * dx;
    p2->vy += dv2 * dy;
    p2->vz += dv2 * dz;
    p2->vr += b/mR;

    colcounter++;

    if (p2->firstcollision && p2->firstcollision != ev)
    {
        removeevent(p2->firstcollision);
    }

    findcollisions(p1);
    findcollisions(p2);
}


/**************************************************
**                 INITEVENTPOOL
** Creates two first events, and sets up
**************************************************/
void initeventpool()
{
    numeventlists = ceil(maxscheduletime / eventlisttime);
    maxscheduletime = numeventlists * eventlisttime;
    if (numeventlists > MAXNUMEVENTLISTS)
    {
        printf("[EDMC-poly] Number of event lists too large: increase MAXNUMEVENTLISTS to at least %d\n", numeventlists);
        exit(3);
    }
    printf("[EDMC-poly] number of lists: %d\n", numeventlists);


    int i;
    event* e;
    for (i = 0; i < MAXEVENTS; i++)		        //Clear eventpool
    {
        e = &(eventlist[MAXEVENTS - i - 1]);	//Fill in backwards, so the first few events are 'special'
        eventpool[i] = e;					    //This includes root, the write events, in that order
        eventpool[i]->child1 = NULL;			//...  Not really used for now, but it might be useful at some point
        eventpool[i]->child2 = NULL;			//Clear children
        nempty++;						        //All events empty so far
    }
    root = eventpool[--nempty];				    //Create root event
    root->time = -99999999999.99;				//Root event is empty, but makes sure every other event has a parent
    root->type = 200;					        //This makes sure we don't have to keep checking this when adding/removing events
    root->parent = NULL;
    event* writeevent = eventpool[--nempty];	//Pick first unused event
    writeevent->time = 0;
    writeevent->type = 100;
    root->child2 = writeevent;
    writeevent->parent = root;
    printf("[EDMC-poly] Event tree initialized: %d events\n", MAXEVENTS - nempty);
}


/**************************************************
**                  ADDEVENTTOTREE
**************************************************/
void addeventtotree(event* newevent)
{
    double time = newevent->time;
    event* loc = root;
    int busy = 1;
    while (busy)    //Find location to add event into tree (loc)
    {
        if (time < loc->time)   //Go left
        {
            if (loc->child1) loc = loc->child1;
            else
            {
                loc->child1 = newevent;
                busy = 0;
            }
        }
        else    //Go right
        {
            if (loc->child2) loc = loc->child2;
            else
            {
                loc->child2 = newevent;
                busy = 0;
            }
        }
    }
    newevent->parent = loc;
}


/**************************************************
**                  ADDEVENT
**************************************************/
void addevent(event* newevent)
{
    double dt = newevent->time - reftime;

    if (dt < eventlisttime) //Put it in the tree
    {
        newevent->queue = currentlist;
        addeventtotree(newevent);
    }
    else
    {
        int list_id;
        if (dt >= numeventlists * eventlisttime) list_id = numeventlists;   //This also handles int overflow when calculating list_id
        else
        {
            list_id = currentlist + dt / eventlisttime;
            if (list_id >= numeventlists)
            {
                list_id -= numeventlists;
            }
        }

        newevent->queue = list_id;
        newevent->nextq = eventlists[list_id];  //Add to linear list
        newevent->prevq = NULL;
        if (newevent->nextq) newevent->nextq->prevq = newevent;
        eventlists[list_id] = newevent;
    }
}


/**************************************************
**                  CREATEEVENT
**************************************************/
event* createevent(double time, particle* p1, particle* p2, int type)
{
    event* newevent = eventpool[--nempty];		//Pick first unused event
    newevent->time = time;
    newevent->p1 = p1;
    newevent->type = type;
    newevent->p2 = p2;
    //printf("newevent->time: %lf / %lf\n", newevent->time, time);
    addevent(newevent);
    return newevent;
}


/**************************************************
**                     ADDNEXTEVENTLIST
**************************************************/
void addnexteventlist()
{
    do
    {
        currentlist++;
        if (currentlist == numeventlists) currentlist = 0;
        reftime += eventlisttime;
    } while (eventlists[currentlist] == NULL);

    //   printf("Currentlist is now %d (%lf)\n", currentlist, reftime);

    event* ev = eventlists[currentlist];
    event* nextev;
    while (ev)
    {
        nextev = ev->nextq;
        //         if (ev->type != 0 || ev->counter2 == ev->p2->counter) 
        addeventtotree(ev);
        ev = nextev;
        listcounter1++;
    }
    eventlists[currentlist] = NULL;
    ev = eventlists[numeventlists];//Overflow queue
    eventlists[numeventlists] = NULL;
    while (ev)
    {
        nextev = ev->nextq;
        addevent(ev);
        ev = nextev;
        listcounter2++;
    }
    mergecounter++;
}


/**************************************************
**                  REMOVEEVENT
**************************************************/
void removeevent(event* oldevent)
{

    //event* ev = oldevent;
    //if (ev->type != 8) printf("Removing event: %lf, ev: %d, part: %d, %d\n", ev->time, ev->type, ev->p1 ? ev->p1->number : -1, ev->p2 ? ev->p2->number : -1);

    if (oldevent->queue != currentlist)
    {
        if (oldevent->nextq) oldevent->nextq->prevq = oldevent->prevq;
        if (oldevent->prevq) oldevent->prevq->nextq = oldevent->nextq;
        else
        {
            eventlists[oldevent->queue] = oldevent->nextq;
        }
        eventpool[nempty++] = oldevent;     //Put the removed event back in the event pool.
        return;
    }

    event* parent = oldevent->parent;
    event* node;					//This node will be attached to parent in the end


    if (oldevent->child1 == NULL)			//Only one child: easy to delete
    {
        node = oldevent->child2;			//Child2 is attached to parent
        if (node)
        {
            node->parent = parent;
            oldevent->child2 = NULL;			//Clear child, so createevent doesn't have to do it
        }
    }
    else if (oldevent->child2 == NULL)		//Only one child again
    {
        node = oldevent->child1;			//Child1 is attached to parent
        node->parent = parent;
        oldevent->child1 = NULL;
    }
    else	//Node to delete has 2 children
    {       //In this case: a) Find first node after oldevent     (This node will have no child1)
            //              b) Remove this node from the tree     (Attach node->child2 to node->parent)
            //              c) Put this node in place of oldevent (Oldevent's children are adopted by node)
        node = oldevent->child2;
        while (node->child1) node = node->child1;	//Find first node of right tree of descendants of oldevent
        event* pnode = node->parent;
        if (pnode != oldevent)			        //node is not a child of oldevent
        {						                //Both of oldevent's children should be connected to node
            pnode->child1 = node->child2;		//Remove node from right tree
            if (node->child2) node->child2->parent = pnode;
            oldevent->child1->parent = node;
            node->child1 = oldevent->child1;
            oldevent->child2->parent = node;
            node->child2 = oldevent->child2;
        }
        else	//This means node == oldevent->child2
        {		//Only child1 has to be attached to node
            oldevent->child1->parent = node;
            node->child1 = oldevent->child1;
        }
        node->parent = parent;
        oldevent->child1 = NULL;
        oldevent->child2 = NULL;
    }
    if (parent->child1 == oldevent) parent->child1 = node;
    else                            parent->child2 = node;
    eventpool[nempty++] = oldevent;     //Put the removed event back in the event pool.
}


/**************************************************
**                  SHOWTREE
** Gives a rough view of the event tree.
** Not so useful except for very small trees
**************************************************/
void showtree()
{
    shownode(root);
}

void shownode(event* ev)
{
    int c1 = 0, c2 = 0, p = 0;
    if (ev->child1) c1 = ev->child1->type;
    if (ev->child2) c2 = ev->child2->type;
    if (ev->parent) p = ev->parent->type;
    printf("%3d => %3d => %d (p: %d, %d)\n", p, ev->type, c1, ev->p1->number, ev->p2->number);
    printf("           => %d \n", c2);

    if (ev->child1) shownode(ev->child1);
    if (ev->child2) shownode(ev->child2);
}


/**************************************************
**                  CHECKTREE
** Checks the tree for possible errors.
**  1 ) node's parent doesn't point to node
**  1b) node->t > parent->t
**  1c) node->t < parent->t
**  2 ) A non-root node lacks a parent
**  3 ) node's child1 doesn't point to parent
**  3b) node->t < child1->t
**  4 ) node's child2 doesn't point to parent
**  4b) node->t > child2->t
** Also checks if all events are in the tree
**************************************************/
void checktree()
{
    return;
    int t = checknode(root);
    if (t != MAXEVENTS - nempty) printf("Error: %d, %d\n", t, MAXEVENTS - nempty);
}

int checknode(event* node)
{
    static int count = 0;
    if (node->parent)
    {
        if (node->parent->child1 != node && node->parent->child2 != node) printf("Error 1\n");
        if (node->parent->child1 == node && node->time > node->parent->time) printf("Error 1b\n");
        if (node->parent->child2 == node && node->time < node->parent->time) printf("Error 1c\n");
    }
    else
    {
        if (root != node) printf("Error 2\n");
        count = 0;
    }
    if (node->child1)
    {
        checknode(node->child1);
        if (node->child1->parent != node) printf("Error 3\n");
        if (node->child1->time > node->time) printf("Error 3b\n");
    }
    if (node->child2)
    {
        checknode(node->child2);
        if (node->child2->parent != node) printf("Error 4\n");
        if (node->child2->time < node->time) printf("Error 4b\n");
    }
    count++;
    return count;
}


/**************************************************
**                    OVERLAP
** Checks for overlaps without any optimization
**************************************************/
int overlap(particle* part)
{
    particle* p;
    double dx, dy, dz, r2, rm;
    int i;
    double dl = pow(10.0, -10);

    for (i = 0; i < N; i++)
    {
        if (i == part->number) continue;
        p = &(particles[i]);
        dx = part->x - p->x;
        dy = part->y - p->y;
        dz = part->z - p->z;
        if (dx > 0.5 * xsize) dx -= xsize; else if (dx < -0.5 * xsize) dx += xsize;  //periodic boundaries
        if (dy > 0.5 * ysize) dy -= ysize; else if (dy < -0.5 * ysize) dy += ysize;
        if (dz > 0.5 * zsize) dz -= zsize; else if (dz < -0.5 * zsize) dz += zsize;
        r2 = dx * dx + dy * dy + dz * dz;
        rm = p->r + part->r;
        if (r2 < rm * rm - dl)
        {
            printf("[EDMC-poly] ERROR: Overlap: %lf, %d, %d\n", r2, part->number, p->number);
            return 1;
        }
    }
    return 0;
}


/**************************************************
**                    OVERLAPLIST
** Checks for overlaps
** if error is one, allow a small margin of error
**************************************************/
int overlaplist(particle* part, int error)
{
    int cdx, cdy, cdz, cellx, celly, cellz, num;
    particle* p;
    double dx, dy, dz, r2, rm;
    double dl = error * pow(10.0, -10);

    cellx = part->cellx + cx;
    celly = part->celly + cy;
    cellz = part->cellz + cz;
    num = part->number;

    for (cdx = cellx - 1; cdx < cellx + 2; cdx++)
        for (cdy = celly - 1; cdy < celly + 2; cdy++)
            for (cdz = cellz - 1; cdz < cellz + 2; cdz++)
            {
                p = celllist[cdx % cx][cdy % cy][cdz % cz];
                while (p)
                {
                    if (p->number != num)
                    {
                        dx = part->x - p->x;
                        dy = part->y - p->y;
                        dz = part->z - p->z;
                        if (dx > 0.5 * xsize) dx -= xsize; else if (dx < -0.5 * xsize) dx += xsize;  //periodic boundaries
                        if (dy > 0.5 * ysize) dy -= ysize; else if (dy < -0.5 * ysize) dy += ysize;
                        if (dz > 0.5 * zsize) dz -= zsize; else if (dz < -0.5 * zsize) dz += zsize;
                        r2 = dx * dx + dy * dy + dz * dz;
                        rm = p->r + part->r;
                        if (r2 < rm * rm - dl)
                        {
                            //            printf ("Overlap: %lf, %d, %d\n", r2, part->number, p->number);
                            return 1;
                        }
                    }
                    p = p->next;
                }
            }
    return 0;
}


/**************************************************
**                    WRITE
** Writes a movie
**************************************************/
void write()
{
    static int counter = 0;
    static int first = 1;
    static double lastsnapshottime = -999999999.9;
    static double lasttime = 0;

    int i;
    particle* p;
    FILE* file;

    double en = 0, enR = 0, packfrac = 0;
    int minneigh = 1000000, maxneigh = 0;
    double minr = maxparticleradius, maxr = minparticleradius;
    for (i = 0; i < N; i++)
    {
        p = particles + i;
        update(p);
        en += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
        enR += p->vr*p->vr;
        packfrac += p->r*p->r*p->r;
        if(p->nneigh < minneigh) minneigh = p->nneigh;
        if(p->nneigh > maxneigh) maxneigh = p->nneigh;
        if(p->r > maxr) maxr = p->r;
        if(p->r < minr) minr = p->r;
    }
    double temperature = 0.5 * en / (double)N / 1.5;
    double totkinen = (en+enR)/ (double)N;
    double volume = xsize * ysize * zsize;
    packfrac *= 4.0/3.0 * M_PI / volume;

    // checktree();
    // checkcells();
    double dens = N / volume;
    double timeint = time - lasttime;
    double press = -dvtot / (3.0 * volume * timeint);
    double pressid = dens;
    double presstot = colcounter ? press + pressid : 0;
    dvtot = 0;

    //Stress tensor
    double Pxx = pressid - Sij[0] / (volume * timeint);
    double Pyy = pressid - Sij[1] / (volume * timeint);
    double Pzz = pressid - Sij[2] / (volume * timeint);
    Sij[0] = 0.0f;
    Sij[1] = 0.0f;
    Sij[2] = 0.0f;

    lasttime = time;

    double listsize1 = (double)listcounter1 / mergecounter;
    double listsize2 = (double)listcounter2 / mergecounter;
    if (mergecounter == 0) listsize1 = listsize2 = 0.0;
    listcounter1 = listcounter2 = mergecounter = 0;

    printf("[EDMC-poly] Simtime: %lf, Collisions: %u, Packfrac: %lf, Press: %lf, Pxx: %lf, Pyy: %lf, Pzz: %lf, T: %lf (%lf), Listsizes: (%lf, %lf), Neighbors: %d - %d\n", 
                time, colcounter, packfrac, presstot, Pxx, Pyy, Pzz, temperature, totkinen, listsize1, listsize2, minneigh, maxneigh);
    char filename[200];
    if (makesnapshots && time - lastsnapshottime > snapshotinterval - 0.001)
    {
        if (usekickback) sprintf(filename, "mov.n%d.v%.4lf.b%.4lf.osph", N, xsize*ysize*zsize, binsize);
        else sprintf(filename, "mov.n%d.v%.4lf.osph", N, xsize * ysize * zsize);
        if (first) { first = 0; file = fopen(filename, "w"); }
        else                     file = fopen(filename, "a");
        fprintf(file, "%d\n%.12lf %.12lf %.12lf\n", (int)N, xsize, ysize, zsize);
        for (i = 0; i < N; i++)
        {
            p = &(particles[i]);

            if (maxr == minr) maxr = minr+0.01;
            double colorpara = (p->r - minr) / (maxr - minr);

            fprintf(file, "%c %.12lf  %.12lf  %.12lf  %lf %lf\n", 'a' + p->type, p->x + xsize * p->boxestraveledx, p->y + ysize * p->boxestraveledy, p->z + zsize * p->boxestraveledz, p->r, colorpara);
        }
        fclose(file);
        lastsnapshottime = time;
    }
    if (verbose && vfile!=NULL && timeint > 0){
        fprintf(vfile, "%lf %lf %lf %lf %lf %lf %lf %lf\n", time, packfrac, presstot, Pxx, Pyy, Pzz, temperature, totkinen);
    }
    if (temperature > 1.5) thermostatinterval *= 0.5;

    if (usekickback && binfile!=NULL){
        fprintf(binfile, "%d\n%.12lf %.12lf %.12lf %.12lf\n", (int)N, xsize, ysize, zsize, binsize);
        for (i = 0; i < N; i++){
            p = &(particles[i]);
            fprintf(binfile, "%d %.6lf %.6lf\n", p->binnumber, p->r, (qbar+i*(2*6+1)+6)->re);
        }
    }

    counter++;
    colcounter = 0;

    createevent(time + writeinterval, NULL, NULL, 100);     //Add next write event
}


/**************************************************
**                    CHECKCELLS
** Checks the cell list for possible errors
**************************************************/
void checkcells()
{
    int x, y, z, count = 0, on = 0, errorcount = 0;
    double dl = 0.0000001;
    particle* p;
    for (x = 0; x < cx; x++)
        for (y = 0; y < cy; y++)
            for (z = 0; z < cz; z++)
            {
                p = celllist[x][y][z];
                if (p)
                {
                    if (p->prev) printf("[EDMC-poly] First part has a prev (%d, %d, %d) %d, %d\n",
                        x, y, z, p->number, p->prev->number);
                    while (p)
                    {
                        on = 0;
                        count++;
                        if (p->cellx != x || p->celly != y || p->cellz != z)
                        {
                            printf("[EDMC-poly] Cell error: %d, %d, %d / %d, %d, %d\n", x, y, z, p->cellx, p->celly, p->cellz);
                            exit(3);
                        }
                        if (p->x < cxsize * x - dl || p->x > cxsize * (x + 1) + dl)
                        {
                            printf("[EDMC-poly] wrong cell x: %lf, %d, %d, %lf\n", p->x, x, p->number, p->vx);
                            printf("[EDMC-poly] wrong cell x: %lf ? %lf ? %lf\n", cxsize*x-dl, p->x, cxsize*(x+1)+dl);
                            on = 1;
                            exit(3);
                        }
                        if (p->y < cysize * y - dl || p->y > cysize * (y + 1) + dl)
                        {
                            printf("[EDMC-poly] wrong cell y: %lf, %d, %d, %lf\n", p->y, y, p->number, p->vy);
                            printf("[EDMC-poly] wrong cell y: %lf ? %lf ? %lf\n", cysize*y-dl, p->y, cxsize*(y+1)+dl);
                            on = 1;
                            exit(3);
                        }
                        if (p->z < czsize * z - dl || p->z > czsize * (z + 1) + dl)
                        {
                            printf("[EDMC-poly] wrong cell z: %lf, %d, %d, %lf\n", p->z, z, p->number, p->vz);
                            printf("[EDMC-poly] wrong cell z: %lf ? %lf ? %lf\n", czsize*z-dl, p->z, czsize*(z+1)+dl);
                            on = 1;
                            exit(3);
                        }
                        if (on) errorcount++;
                        if (p->next)
                        {
                            if (p->next->prev != p) printf("[EDMC-poly] link error: %d, %d, %d\n",
                                p->number, p->next->number, p->next->prev->number);
                        }
                        p = p->next;
                    }
                }
            }
    if (count != N) printf("[EDMC-poly] error in number of particles (%d)\n", count);
    printf("[EDMC-poly] errorcount: %d / %d\n", errorcount, N);
}


/**************************************************
**                    BACKINBOX
** Just for initialization
**************************************************/
void backinbox(particle* p)
{
    p->x -= xsize * floor(p->x / xsize);
    p->y -= ysize * floor(p->y / ysize);
    p->z -= zsize * floor(p->z / zsize);
}


/**************************************************
**                    THERMOSTAT
**************************************************/
void thermostat(event* ev)
{
    printf("[EDMC-poly] therm\n");
    if (ev)
    {
        int i, num;
        particle* p;
        time = ev->time;
        int freq = N / 100;
        if (freq == 0) freq = 1;
        for (i = 0; i < freq; i++)
        {
            num = genrand_real2() * N;			//Random particle
            p = particles + num;
            double imsq = 1.0 / sqrt(p->mass);
            update(p);
            p->vx = random_gaussian() * imsq;			//Kick it
            p->vy = random_gaussian() * imsq;
            p->vz = random_gaussian() * imsq;
            p->counter++;
            removeevent(p->firstcollision);
            findcollisions(p);
        }
        removeevent(ev);
    }
    createevent(time + thermostatinterval, NULL, NULL, 200);     //Add next write interval
}


/**************************************************
**                    WRITELAST
**************************************************/
void writelast()
{
    int i;
    particle* p;
    FILE* file;
    char filename[200];
    sprintf(filename, "last.sph");
    file = fopen(filename, "w");
    fprintf(file, "%d\n%lf %lf %lf\n", (int)N, xsize, ysize, zsize);
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        fprintf(file, "%c %.12lf  %.12lf  %.12lf  %.12lf\n", 'a' + p->type, p->x, p->y, p->z, p->r);
    }
    fclose(file);
}


/**************************************************
**                  RANDOM_GAUSSIAN 
**************************************************/
double random_gaussian()
{
    static int have_deviate = 0;
    static double u1, u2;
    double  x1, x2, w;

    if (have_deviate)
    {
        have_deviate = 0;
        return u2;
    }
    else
    {
        do
        {
            x1 = 2.0 * genrand_real2() - 1.0;
            x2 = 2.0 * genrand_real2() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);
        w = sqrt((-2.0 * log(w)) / w);
        u1 = x1 * w;
        u2 = x2 * w;
        have_deviate = 1;
        return u1;
    }
}


/**************************************************
**              FIND_LOCMIN_P3
** Identifies the local minimum for 3rd order
** polynomial
**************************************************/
double find_locmin_p3() {
    double delta = 4.0f*a2*a2 - 12.0f*a1*a3;
    if (delta == 0.0f){
        double r1 = (-2.0f*a3) / (6.0f*a3);
        double muppr1 = 2.0f*a2 + 6.0f*a3*r1;
        if (muppr1 > 0.0f){
            printf("[EDMC-poly] No local minimum - inflexion point at rmin: %.4lf - mu is monotonically increasing\n", r1);
            r1 = minparticleradius;
        }
        else{
            printf("[EDMC-poly] No local minimum - inflexion point at rmin: %.4lf - mu is monotonically decreasing\n", r1);
            r1 = maxparticleradius;
        }
        return r1;
    }
    if (delta < 0.0f){
        //Derivative has negative discriminant, presence of an inflexion point
        double monotony = a1 * (minparticleradius - maxparticleradius) +
            a2 * (minparticleradius*minparticleradius - maxparticleradius*maxparticleradius) +
            a3 * (minparticleradius*minparticleradius*minparticleradius - maxparticleradius*maxparticleradius*maxparticleradius);
        double r1 = - a2 / (3.0f * a3);
        if (monotony > 0.0f){
            printf("[EDMC-poly] WARNING: No local minimum - inflexion point at rmin: %.4lf - mu is monotonously decreasing\n", r1);
            r1 = maxparticleradius;
        }
        else{
            printf("[EDMC-poly] WARNING: No local minimum - inflexion point at rmin: %.4lf - mu is monotonously increasing\n", r1);
            r1 = minparticleradius;
        }
        return r1;
    }
    double r1 = (-2.0f*a2 + sqrt(delta)) / (6.0f*a3);
    double r2 = (-2.0f*a2 - sqrt(delta)) / (6.0f*a3);
    double muppr1 = 2.0f*a2 + 6.0f*a3*r1;
    if (muppr1 > 0.0f){
        if (r1 < minparticleradius){
            r1 = minparticleradius;
            printf("[EDMC-poly] Local minimum found outside of considered particle size range: rmin: %.4lf\n", r1);
        }
        if (r1 > maxparticleradius){
            r1 = maxparticleradius;
            printf("[EDMC-poly] Local minimum found outside of considered particle size range: rmin: %.4lf\n", r1);
        }
        return r1;
    }
    else{
        if (r2 < minparticleradius){
            r2 = minparticleradius;
            printf("[EDMC-poly] Local minimum found outside of considered particle size range: rmax: %.4lf\n", r2);
        }
        if (r2 > maxparticleradius){
            r2 = maxparticleradius;
            printf("[EDMC-poly] Local minimum found outside of considered particle size range: rmax: %.4lf\n", r2);
        }
        return r2;
    }
}


/**************************************************
**              FIND_REALROOTS_P3
** Finds the real roots of a third order polynomial
** with real coefficients using the trigonometric
** solution
**************************************************/
double find_realroots_p3(double mu, int k){
    double depressed_q = (2.0f*a2*a2*a2 - 9.0f*a1*a2*a3 + 27.0f*a3*a3*(a0-mu)) / (27.0f*a3*a3*a3);
    double argarccos = (3.0f*depressed_q) / (2.0f*depressed_p) * sqrt(-3.0f/depressed_p);
    double argcos = acos(argarccos) / 3.0f - (2.0f*M_PI*(double)k/3.0f);
    double t = 2.0f*sqrt(-depressed_p / 3.0f) * cos(argcos);
    return t - a2 / (3.0f*a3);
}


/**************************************************
**          FINDRADIALCOLLISION_TAILORED
** Tailored field interaction
**************************************************/
double findradialcollision_tailored(particle* p1, double tmin)
{
    double rcur = p1->r;
    double vrcur = p1->vr;
    double rnew = rcur;
    double rtemp;
    int nroots = 0;
    double factor1 = a2*a2 - 3.0f*a3*a1;

    if(vrcur > 0){
        //Currently growing
        double ucurrent = umin;
        if(rcur > rmin) ucurrent = a0 + a1*rcur + a2*rcur*rcur + a3*rcur*rcur*rcur;
        else if (rmin == maxparticleradius){
            //Specific case of monotously decreasing chem. pot.
            double t = (maxparticleradius - rcur) / vrcur;
            if (t<0.0f) printf("[EDMC-poly] WARNING: t<0 (%.12lf) (gmd) - ucur: %.6lf rcur: %.6lf unew: %.6lf rnew: %.6lf\n",t , ucurrent, rcur, umin, maxparticleradius);
            return t;
        }
        double unew = ucurrent + random_exponential();      //Energy at collision
        //Identify how many roots
        double factor2 = 2.0f*a2*a2*a2 - 9.0f*a1*a2*a3 + 27.0f*a3*a3*(a0-unew);
        double Delta = (4.0f*factor1*factor1*factor1 - factor2*factor2) / (27.0f*a3*a3);
        if (Delta > 0.0f){
            //Three real roots - trigonometric solution
            for (int i=0; i<3; i++){
                rtemp = find_realroots_p3(unew, i);
                if (rtemp > rmin){
                    //Sliding down is free: only compare w.r.t. rmin
                    roots[nroots] = rtemp;
                    nroots++;
                }
            }
            if (nroots == 0){
                printf("[EDMC-poly] WARNING: no suitable root found (g, D>0)\n");
                rnew = maxparticleradius;
            }
            else{
                rnew = roots[0];
                for (int i=1; i<nroots; i++){
                    if (roots[i] < rnew) rnew = roots[i];
                }
            }
        }
        else if (Delta == 0.0f){
            double delta = a2*a2 - 3.0f*a3*a1; 
            if (delta == 0){
                //Triple root
                rnew = -a2/(3.0f*a3);
            }
            else{
                //Double root + simple root
                rnew = (9.0f*(a0-unew)*a3 - a1*a2) / (2.0f*(a2*a2 - 3.0f*a1*a3));
                rtemp = (4.0f*a1*a2*a3 - 9.0f*a3*a3*(a0-unew) - a2*a2*a2) / (a3*(a2*a2 - 3.0f*a1*a3));
                if (rnew < rcur){
                    if (rtemp < rcur){
                        rnew = rcur;
                        printf("[EDMC-poly] WARNING: no suitable root found (g, D=0)\n");
                    }
                    else rnew = rtemp;
                }
                else{
                    if (rtemp > rcur && rtemp < rnew) rnew = rtemp;
                }
            }
        }
        else{
            //One real root - hyperbolic solution
            double depressed_q = (2.0f*a2*a2*a2 - 9.0f*a1*a2*a3 + 27.0f*a3*a3*(a0-unew)) / (27.0f*a3*a3*a3);
            double delta = 4.0f*depressed_p*depressed_p*depressed_p + 27.0f*depressed_q*depressed_q;
            double sign_q = 1.0f;
            if (delta > 0.0f && depressed_p < 0.0f){
                if (depressed_q < 0.0f) sign_q = -1.0f;
                double argarccosh = (-3.0f*sign_q*depressed_q) / (2.0f*depressed_p) * sqrt(-3.0f/depressed_p);
                double argcosh = acosh(argarccosh) / 3.0f;
                rnew = -2.0f * sign_q * sqrt(-depressed_p/3.0f) * cosh(argcosh) - a2/(3.0f*a3);
            }
            else if (depressed_p > 0){
                double argarcsinh = (3.0f*depressed_q) / (2.0f*depressed_p) * sqrt(3.0f/depressed_p);
                double argsinh = asinh(argarcsinh) / 3.0f;
                rnew = -2.0f * sqrt(depressed_p / 3.0f) * sinh(argsinh) - a2/(3.0f*a3);
            }
            else{
                rnew = rcur;
                printf("[EDMC-poly] WARNING: no suitable root found (g, D<0)\n");
            }
        }
        if (rnew > maxparticleradius) rnew = maxparticleradius;
        double t = (rnew - rcur) / vrcur;
        if (t<0.0f) printf("[EDMC-poly] WARNING: t<0 (%.12lf) (g) - ucur: %.6lf rcur: %.6lf unew: %.6lf rnew: %.6lf\n",t , ucurrent, rcur, unew, rnew);
        return t;
    }
    else{
        //Currently shrinking
        double ucurrent = umin;
        if(rcur < rmin) ucurrent = a0 + a1*rcur + a2*rcur*rcur + a3*rcur*rcur*rcur;
        else if (rmin == minparticleradius){
            //Specific case of monotously increasing chem. pot.
            double t = (minparticleradius - rcur) / vrcur;
            if (t<0.0f) printf("[EDMC-ply] WARNING: t<0 (%.12lf) (gmi) - ucur: %.6lf rcur: %.6lf unew: %.6lf rnew: %.6lf\n",t , ucurrent, rcur, umin, minparticleradius);
            return t;
        }
        double unew = ucurrent + random_exponential();      //Energy at collision
        //Identify how many roots
        double factor2 = 2.0f*a2*a2*a2 - 9.0f*a1*a2*a3 + 27.0f*a3*a3*(a0-unew);
        double Delta = (4.0f*factor1*factor1*factor1 - factor2*factor2) / (27.0f*a3*a3);
        if (Delta > 0.0f){
            //Three real roots - trigonometric solution
            for (int i=0; i<3; i++){
                rtemp = find_realroots_p3(unew, i);
                if (rtemp < rmin){
                    //Sliding down is free: only compare w.r.t. rmin
                    roots[nroots] = rtemp;
                    nroots++;
                }
            }
            if (nroots == 0){
                printf("[EDMC-poly] WARNING: no suitable root found (s, D>0)\n");
                rnew = minparticleradius;
            }
            else{
                rnew = roots[0];
                for (int i=1; i<nroots; i++){
                    if (roots[i] > rnew) rnew = roots[i];
                }
            }
        }
        else if (Delta == 0.0f){
            double delta = a2*a2 - 3.0f*a3*a1;
            if (delta == 0.0f){
                //Triple root
                rnew = -a2/(3.0f*a3);
            }
            else{
                //Double root + simple root
                rnew = (9.0f*(a0-unew)*a3 - a1*a2) / (2.0f*(a2*a2 - 3.0f*a1*a3));
                rtemp = (4.0f*a1*a2*a3 - 9.0f*a3*a3*(a0-unew) - a2*a2*a2) / (a3*(a2*a2 - 3.0f*a1*a3));
                if (rnew > rcur){
                    if (rtemp > rcur){
                        rnew = rcur;
                        printf("[EDMC-poly] WARNING: no suitable root found (s, D=0)\n");
                    }
                    else rnew = rtemp;
                }
                else if (rtemp < rcur && rtemp > rnew) rnew = rtemp;
            }
        }
        else{
            //One real root - hyperbolic solution
            double depressed_q = (2.0f*a2*a2*a2 - 9.0f*a1*a2*a3 + 27.0f*a3*a3*(a0-unew)) / (27.0f*a3*a3*a3);
            double delta = 4.0f*depressed_p*depressed_p*depressed_p + 27.0f*depressed_q*depressed_q;
            double sign_q = 1.0f;
            if (delta > 0.0f && depressed_p < 0.0f){
                if (depressed_q < 0.0f) sign_q = -1.0f;
                double argarccosh = (-3.0f*sign_q*depressed_q) / (2.0f*depressed_p) * sqrt(-3.0f/depressed_p);
                double argcosh = acosh(argarccosh) / 3.0f;
                rnew = -2.0f * sign_q * sqrt(-depressed_p/3.0f) * cosh(argcosh) - a2/(3.0f*a3);
            }
            else if (depressed_p > 0){
                double argarcsinh = (3.0f*depressed_q) / (2.0f*depressed_p) * sqrt(3.0f/depressed_p);
                double argsinh = asinh(argarcsinh) / 3.0f;
                rnew = -2.0f * sqrt(depressed_p / 3.0f) * sinh(argsinh) - a2/(3.0f*a3);
            }
            else{
                rnew = rcur;
                printf("[EDMC-poly] WARNING: no suitable root found (s, D<0)\n");
            }
        }
        if (rnew < minparticleradius) rnew = minparticleradius;
        double t = (rnew - rcur) / vrcur;
        if (t < 0.0f) printf("[EDMC-poly] WARNING: t<0 (s) - ucur: %.6lf rcur: %.6lf unew: %.6lf rnew: %.6lf\n", ucurrent, rcur, unew, rnew);
        return t;        
    }
}




/**************************************************
**                    FCCPARTICLES
** Positions particles on an FCC latice in a cubic
** box
** Particles start small
** NOTE: only tolerates N such that cbrt(N) is an
**       integer
**************************************************/
void fccparticles(){
    double vol = N / density;
    printf("[EDMC-poly] Volume: %lf\n", vol);

    xsize = cbrt(vol);
    ysize = xsize;
    zsize = ysize;
    initcelllist();
    particle* p;
    
    int cbrtN = cbrt(N);
    int n = 0;
    if ((double)cbrtN / cbrt(4.0f) * sqrt(2.0f) > xsize){
        printf("[EDMC-poly] ERROR at FCC initialization (will be fixed...)\n");
        exit(3);
    }
    double interspace = cbrt(4.0f / density);
    for (int i=0; i<cbrtN/cbrt(4.0f); i++){     //First put particles on a fcc lattice
        for (int j=0; j<cbrtN/cbrt(4.0f); j++){
            for (int k=0; k<cbrtN/cbrt(4.0f); k++){
                p = &(particles[n]);
                p->x = k*interspace;
                p->y = j*interspace;
                p->z = i*interspace;            // {0,0,0}
                p = &(particles[n+1]);
                p->x = (k+0.5f)*interspace;
                p->y = (j+0.5f)*interspace;
                p->z = i*interspace;            // {a/2, a/2, 0}
                p = &(particles[n+2]);
                p->x = (k+0.5f)*interspace;
                p->y = j*interspace;
                p->z = (i+0.5f)*interspace;     // {a/2, 0, a/2}
                p = &(particles[n+3]);
                p->x = k*interspace;
                p->y = (j+0.5f)*interspace;
                p->z = (i+0.5f)*interspace;     // {0, a/2, a/2}
                n+=4;
            }
        }
    }				

    for (int i=0; i<N; i++){
        p = &(particles[i]);
        p->rtarget = 0.5;
        p->r = 0.8 * p->rtarget;                //Start particles off small, so it's easy to make a random configuration
        p->mass = 1;
        p->type = 0;
        p->number = i;
        do{
            p->cellx = p->x / cxsize;			//Find particle's cell
            p->celly = p->y / cysize;
            p->cellz = p->z / czsize;
        } while (overlaplist(p, 0));
        double sqm = 1.0 / sqrt(p->mass);
        p->xn = p->x; p->yn = p->y; p->zn = p->z;
        p->vx = (genrand_real2() - 0.5) / sqm;
        p->vy = (genrand_real2() - 0.5) / sqm;
        p->vz = (genrand_real2() - 0.5) / sqm;
        p->t = 0;   //r and v known at t=0
        p->next = celllist[p->cellx][p->celly][p->cellz];	//Add particle to celllist
        if (p->next) p->next->prev = p;			//Link up list
        celllist[p->cellx][p->celly][p->cellz] = p;
        p->prev = NULL;
    }
}


/**************************************************
**              FCCAPARTICLES_COEX
** Positions particles on an FCC latice in an
** elongated box for use in direct coexistence 
** semigrand simulations
** Crystal orientation follows the FCCa convention
** Particles start small, fluid particles are
** inserted to amount up to N particles
**************************************************/
void fccaparticles_coex(){
    if (density == 0.0f) density = (fdensity + Xdensity) / 2.0f;
    int n = 1.0f/FCClatticeparam * cbrt((double)N/(3.0f*density));
    xsize = n*FCClatticeparam;
    ysize = xsize;
    zsize = (double)N / (density*xsize*xsize);
    int m = propVolX * zsize / FCClatticeparam;
    double zsizeXinit = m*FCClatticeparam;
    initcelllist();

    void setup_part(int partid){
        particle* p = &(particles[partid]);
        p->cellx = p->x / cxsize;
        p->celly = p->y / cysize;
        p->cellz = p->z / czsize;
        p->rtarget = 0.5;
        p->r = 0.9 * p->rtarget; 
        p->mass = 1;
        p->type = 0;
        p->number = partid;
        double sqm = 1.0 / sqrt(p->mass);
        p->xn = p->x; p->yn = p->y; p->zn = p->z;
        p->vx = (genrand_real2() - 0.5) / sqm;
        p->vy = (genrand_real2() - 0.5) / sqm;
        p->vz = (genrand_real2() - 0.5) / sqm;
        p->t = 0;						//r and v known at t=0
        p->next = celllist[p->cellx][p->celly][p->cellz];	//Add particle to celllist
        if (p->next) p->next->prev = p;			//Link up list
        celllist[p->cellx][p->celly][p->cellz] = p;
        p->prev = NULL;
    }

    // Build crystal lattice
    particle* p;
    int pid = 0;
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                p = &(particles[pid]);
                p->x = k*FCClatticeparam;
                p->y = j*FCClatticeparam;
                p->z = i*FCClatticeparam;            // {0,0,0}
                setup_part(pid);
                p = &(particles[pid+1]);
                p->x = (k+0.5f)*FCClatticeparam;
                p->y = (j+0.5f)*FCClatticeparam;
                p->z = i*FCClatticeparam;            // {a/2, a/2, 0}
                setup_part(pid+1);
                p = &(particles[pid+2]);
                p->x = (k+0.5f)*FCClatticeparam;
                p->y = j*FCClatticeparam;
                p->z = (i+0.5f)*FCClatticeparam;     // {a/2, 0, a/2}
                setup_part(pid+2);
                p = &(particles[pid+3]);
                p->x = k*FCClatticeparam;
                p->y = (j+0.5f)*FCClatticeparam;
                p->z = (i+0.5f)*FCClatticeparam;     // {0, a/2, a/2}
                setup_part(pid+3);
                pid+=4;
            }
        }
    }
   
    // Insert remaining particles
    for (int i=pid; i<N; i++){
        p = &(particles[i]);
        p->type = 0;
        p->rtarget = 0.5f;
        p->r = 0.6 * p->rtarget;    //Start particles off small, so it's easy to make a random configuration
        p->mass = 1.0f;
        p->type = 0;
        p->number = i;
        do{
            p->x = genrand_real2() * xsize;			//Random location
            p->y = genrand_real2() * ysize;
            p->z = genrand_real2() * zsize;
            p->cellx = p->x / cxsize;
            p->celly = p->y / cysize;
            p->cellz = p->z / czsize;
        } while (overlaplist(p, 0) || p->z <= zsizeXinit);
        double sqm = 1.0 / sqrt(p->mass);
        p->xn = p->x; p->yn = p->y; p->zn = p->z;
        p->vx = (genrand_real2() - 0.5) / sqm;
        p->vy = (genrand_real2() - 0.5) / sqm;
        p->vz = (genrand_real2() - 0.5) / sqm;
        p->t = 0;						//r and v known at t=0
        p->next = celllist[p->cellx][p->celly][p->cellz];	//Add particle to celllist
        if (p->next) p->next->prev = p;			//Link up list
        celllist[p->cellx][p->celly][p->cellz] = p;
        p->prev = NULL;
    }
}


/**************************************************
**              FCCBPARTICLES_COEX
** Positions particles on an FCC latice in an
** elongated box for use in direct coexistence 
** semigrand simulations
** Crystal orientation follows the FCCb convention
** Particles start small, fluid particles are
** inserted to amount up to N particles
**************************************************/
void fccbparticles_coex(){
    if (density == 0.0f) density = (fdensity + Xdensity) / 2.0f;
    double xsizetemp = cbrt((double)N/(3.0f*density));
    int m = sqrt(2.0f)*xsizetemp / (FCClatticeparam*sqrt(3.0f));
    xsize = (double)m*FCClatticeparam*sqrt(3.0f)/sqrt(2.0f);
    int n = xsizetemp / (FCClatticeparam*sqrt(3.0f));
    ysize = (double)n*FCClatticeparam*sqrt(3.0f);
    zsize = (double)N / (density*xsize*ysize);
    int o = propVolX * zsize * sqrt(2.0f) / FCClatticeparam;
    double zsizeXinit = o*FCClatticeparam/sqrt(2.0f);
    double ax = FCClatticeparam*sqrt(3.0f)/sqrt(2.0f);
    double ay = FCClatticeparam*sqrt(3.0f);
    double az = FCClatticeparam/sqrt(2.0f);
    initcelllist();

    void setup_part(int partid){
        particle* p = &(particles[partid]);
        p->cellx = p->x / cxsize;
        p->celly = p->y / cysize;
        p->cellz = p->z / czsize;
        p->rtarget = 0.5;
        p->r = 0.9 * p->rtarget; 
        p->mass = 1;
        p->type = 0;
        p->number = partid;
        double sqm = 1.0 / sqrt(p->mass);
        p->xn = p->x; p->yn = p->y; p->zn = p->z;
        p->vx = (genrand_real2() - 0.5) / sqm;
        p->vy = (genrand_real2() - 0.5) / sqm;
        p->vz = (genrand_real2() - 0.5) / sqm;
        p->t = 0;						//r and v known at t=0
        p->next = celllist[p->cellx][p->celly][p->cellz];	//Add particle to celllist
        if (p->next) p->next->prev = p;			//Link up list
        celllist[p->cellx][p->celly][p->cellz] = p;
        p->prev = NULL;
    }

    // Build crystal lattice
    particle* p;
    int pid = 0;
    for (int i=0; i<o; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<m; k++){
                p = &(particles[pid]);      // {0, 0, 0}
                p->x = k*ax;
                p->y = j*ay;
                p->z = i*az;
                setup_part(pid);
                p = &(particles[pid+1]);    // {a*sqrt(3)/(2*sqrt(2)), 0, a/(2*sqrt(2))}
                p->x = (k+0.5f)*ax;
                p->y = j*ay;
                p->z = (i+0.5f)*az;
                setup_part(pid+1);
                p = &(particles[pid+2]);    // {a*sqrt(3)/(6*sqrt(2)), a/sqrt(3), a/(2*sqrt(2))}
                p->x = (k+1.0f/6.0f)*ax;
                p->y = (j+1.0f/3.0f)*ay;
                p->z = (i+0.5f)*az;
                setup_part(pid+2);
                p = &(particles[pid+3]);    // {a*sqrt(2)/sqrt(3), a/sqrt(3), 0}
                p->x = (k+2.0f/3.0f)*ax;
                p->y = (j+1.0f/3.0f)*ay;
                p->z = i*az;
                setup_part(pid+3);
                p = &(particles[pid+4]);    // {a/sqrt(6), 2*a/sqrt(3), 0}
                p->x = (k+1.0f/3.0f)*ax;
                p->y = (j+2.0f/3.0f)*ay;
                p->z = i*az;
                setup_part(pid+4);
                p = &(particles[pid+5]);   // {5*a/(2*sqrt(6)), 2*a/sqrt(3), a/(2*sqrt(2))}
                p->x = (k+5.0f/6.0f)*ax;
                p->y = (j+2.0f/3.0f)*ay;
                p->z = (i+0.5f)*az;
                setup_part(pid+5);
                pid+=6;
            }
        }
    }
   
    // Insert remaining particles
    for (int i=pid; i<N; i++){
        p = &(particles[i]);
        p->type = 0;
        p->rtarget = 0.5f;
        p->r = 0.6 * p->rtarget;    //Start particles off small, so it's easy to make a random configuration
        p->mass = 1.0f;
        p->type = 0;
        p->number = i;
        do{
            p->x = genrand_real2() * xsize;			//Random location
            p->y = genrand_real2() * ysize;
            p->z = genrand_real2() * zsize;
            p->cellx = p->x / cxsize;
            p->celly = p->y / cysize;
            p->cellz = p->z / czsize;
        } while (overlaplist(p, 0) || p->z <= zsizeXinit);
        double sqm = 1.0 / sqrt(p->mass);
        p->xn = p->x; p->yn = p->y; p->zn = p->z;
        p->vx = (genrand_real2() - 0.5) / sqm;
        p->vy = (genrand_real2() - 0.5) / sqm;
        p->vz = (genrand_real2() - 0.5) / sqm;
        p->t = 0;						//r and v known at t=0
        p->next = celllist[p->cellx][p->celly][p->cellz];	//Add particle to celllist
        if (p->next) p->next->prev = p;			//Link up list
        celllist[p->cellx][p->celly][p->cellz] = p;
        p->prev = NULL;
    }
}


/**************************************************
**                  LOADPARTICLES
** Loads initial configuration from a user-supplied
** file
**************************************************/
void loadparticles() {
    char tmp;
    int i, npart;
    particle* p;
    char buffer[255];

    FILE* file;
    file = fopen(inputfilename, "r");
    if (!file) {
        printf("[EDMC-poly] File not found: %s\n", inputfilename);
        exit(3);
    }
    mygetline(buffer, file);
    int ftmp = sscanf(buffer, "%d", &npart);
    if (ftmp != 1) {printf("[EDMC-poly] Read error (N)\n"); exit(3);}
    if (npart != N) {printf("[EDMC-poly] ERROR: number of particles from loaded configuration does not match\n"); exit(3);}
    mygetline(buffer, file);
    ftmp = sscanf(buffer, "%lf %lf %lf\n", &xsize, &ysize, &zsize);
    if (ftmp != 3) {printf("[EDMC-poly] Read error (box size)\n"); exit(3);}
    initcelllist();


    for (i = 0; i < N; i++) {
        p = &(particles[i]);
        mygetline(buffer, file);
        ftmp = sscanf(buffer, "%c %lf  %lf  %lf %lf %*f\n", &tmp, &(p->x), &(p->y), &(p->z), &(p->r));
        backinbox(p);
        if (ftmp != 5) {printf("[EDMC-poly] Read error (particle %d)\n String: %s\n", ftmp, buffer); exit(3); }
        p->rtarget = p->r;
        p->type = tmp - 'a';
        p->mass = 1;
        p->number = i;
        p->cellx = p->x / cxsize;
        p->celly = p->y / cysize;
        p->cellz = p->z / czsize;
        double sqm = 1.0 / sqrt(p->mass);
        p->xn = p->x;
        p->yn = p->y;
        p->zn = p->z;
        p->vx = (genrand_real2() - 0.5) / sqm;
        p->vy = (genrand_real2() - 0.5) / sqm;
        p->vz = (genrand_real2() - 0.5) / sqm;
        p->t = 0;						//r and v known at t=0
        p->next = celllist[p->cellx][p->celly][p->cellz];	//Add particle to celllist
        if (p->next) p->next->prev = p;			//Link up list
        celllist[p->cellx][p->celly][p->cellz] = p;
        p->prev = NULL;
    }
    fclose(file);
}


/**************************************************
**                      KICKBACK
** Kicks the particles back in place according to
** the mean displacement of the crystal bulk, used
** to `pin` coexistence interface in place
** NOTE: Requires for solid-like particles to be
** labelled with something else than `a` type which
** is used for fluid-like
**************************************************/
void kickback(event* ev){
    if (ev){
        static int kickcount = 0;
        int Xpartcount = 0;
        double deltaz;
        time = ev->time;

        kickcount++;
        if (kickcount % (int) (writeinterval / kickinterval) == (int) (writeinterval / kickinterval) - 1){
            buildsann();
            calc_order();
            calc_conn();
        }

        // Calculate mean X displacement
        for (int i=0; i<N; i++){
            particle* p = particles+i;
            if (p->type != 'a'){
                Xpartcount++;
                particle* peq = particleseq+i;
                Xoffset += p->z - peq->z;
            }
        }
        Xoffset /= (double) Xpartcount;
        
        // reattribute particles to their bin
        for (int i=0; i<N; i++){
            particle* p = particles+i;
            deltaz = p->z - Xoffset;
            if (deltaz >= zsize) deltaz -= zsize;
            else if (deltaz < 0) deltaz += zsize;
            int binnumber = (int) (deltaz / binsize);
            if (binnumber < 0) printf("[EDMC-poly] WARNING: binnumber: %d z: %lf offset: %lf\n", binnumber, p->z, Xoffset);
            p->binnumber = binnumber;
        }
    }

    createevent(time + kickinterval, NULL, NULL, 150);
}


/**************************************************
**                  DOTPRODSUM
** Computes the sum of dotproducts for a l-defined
** range of vectors couples
**************************************************/
double dotprodsum(cmplx *vec1, cmplx *vec2, int l){
    double res = 0.0f;
    for(int m = -l; m <= l; m++)
        res += (*(vec1 + m + l)).re * (*(vec2 + m + l)).re
            + (*(vec1 + m + l)).im * (*(vec2 + m + l)).im;
    return res;
}


/**************************************************
**                  MINPOW
** Returns 1.0f if m is even, -1.0f if not
 **************************************************/
double minpow(int m){
    if((m & 1) == 1) return -1.0;
    else return 1.0;
}


/**************************************************
**                  GAMMLN
** Computes factorials using the gamma function
**************************************************/
double gammln(double xx){
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146,
                            -86.50532032941677,
                            24.01409824083091,
                            -1.231739572450155,
                            0.1208650973866179e-2,
                            -0.5395239384953e-5
                            };
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}


/**************************************************
**                      FACS
** Computes (l-m)!/(l+m)!
**************************************************/
double facs(int l, int m){
    static double* fac_table = NULL;
    int a, b;
    if(fac_table == NULL){
        fac_table = malloc((2*l+1) * sizeof(*fac_table));
        for(a = 0; a < 2*l+1; a++){
            b = a - l;
            fac_table[a]= exp(gammln(l - b + 1) - gammln(l + b + 1));
        }
    }
    return fac_table[m+l];
}


/**************************************************
**                  PLGNDR
** Calculates the Legendre function P_{l,m}(x) = 
** (1-x**2)**{m/2} (\frac{d}{dx})**m P_l(x) where 
** P_l(x) is the Legendre polynomial defined for x
** on [-1;1]
**************************************************/
float plgndr(int l, int m, double x){
    double fact,
       pll = 0.0f,
       pmm,
       pmmp1,
       somx2;
    int i, ll;

    //Check for normal computation of Legendre polynoms
    if (m < 0 || m > l || fabs(x) > 1.0)
        printf("[EDMC-poly] WARNING: Bad arguments in routine plgndr %i %i %f\n", l, m, fabs(x));
    pmm = 1.0;
    if (m > 0){
        somx2 = sqrt((1.0 - x) * (1.0 + x));
        fact = 1.0;
        for (i = 1; i <= m; i++){
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m){return pmm;}
    else{
        pmmp1 = x * (2 * m + 1) * pmm;
        if (l == (m + 1)){return pmmp1;}
        else{
            for (ll = m + 2; ll <= l; ll++){
                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}


/**************************************************
**              CALC_HARMONICS
** Calculates spherical harmonics (normalized by
** the size of the particle's neighbor list)
**************************************************/
void calc_harmonics(int l, int n1, int n2, cmplx *res1, cmplx *res2){
    double fc, p, f, s, r, sp, spp, c, cp, cpp;
    int m = 0;
    particle* p1 = particles+n1;
    particle* p2 = particles+n2;

    double nnneigh1 = p1->nnneigh;
    double nnneigh2 = p2->nnneigh;

    //Calculate distances & angles
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    double dz = p1->z - p2->z;
    if (dx > 0.5 * xsize) dx -= xsize; else if (dx < -0.5 * xsize) dx += xsize;  //periodic boundaries
    if (dy > 0.5 * ysize) dy -= ysize; else if (dy < -0.5 * ysize) dy += ysize;
    if (dz > 0.5 * zsize) dz -= zsize; else if (dz < -0.5 * zsize) dz += zsize;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    double dxy = 1.0 / sqrt(dx*dx + dy*dy);
    double z = dz / dist;
    double co = dx * dxy;
    double si = dy * dxy;

    //Computes the spherical harmonics for m = 0
    p = plgndr(l,0,z);
    fc = facs(l,0);
    f = sqrt((2*l+1) * INVFPI * fc);
    r = p*f;
    (res1+0)->re += r / nnneigh1;
    (res1+0)->im += 0;
    (res2+0)->re += r * minpow(l) / nnneigh2; // minpow(6)=1
    (res2+0)->im += 0;

    s=0;
    sp=0;
    c=1;
    cp=0;

    for(m = 1; m <= l; m++){
            //For m > 0
            p = plgndr(l,m,z);
            fc = facs(l,m);
            f = sqrt((2 * l + 1) * INVFPI * fc);
            r = p * f;
            //Chebyshev recursive method for computing cosine of multiple angles 
            cpp = cp;
            cp = c;
            if(m == 1){c = co;}
            else{c = 2.0 * co * cp - cpp;}
            //Chebyshev recursive method for computing sine of multiple angles 
            spp = sp;
            sp = s;
            if(m == 1){s = si;}
            else{s = 2.0 * co * sp - spp;}

            (res1+m)->re += r*c / nnneigh1;
            (res1+m)->im += r*s / nnneigh1;
            (res2+m)->re += r*c / nnneigh2;
            (res2+m)->im += r*s / nnneigh2;

            //For m < 0
            r *= minpow(m);
            (res1-m)->re += r*c / nnneigh1;
            (res1-m)->im += -r*s / nnneigh1;
            (res2-m)->re += r*c / nnneigh2;
            (res2-m)->im += -r*s / nnneigh2;
    }
}


/**************************************************
**                  CALC_ORDER
** ---Calculates Lechner & Dellago averaged BOP---
** Calculates ten Wolde \tilde{q}_{lm}
**************************************************/
void calc_order(){
    cmplx* q1;
    cmplx* q2;
    particle* p1;
    particle* p2;
    double temp;
    const int l = 6;
    memset(bop6, (int) 0.0, sizeof(*bop6) * N * (l * 2 + 1));
    memset(qbar, (int) 0.0, sizeof(*qbar) * N * (l * 2 + 1));

    //Compute q_{lm}(i)
    for(int i = 0; i < N; i++){
        p1 = particles + i;
        q1 = (bop6 + i * (2 * l + 1) + l);
        for(int j = 0; j < p1->nnneigh; j++){
            p2 = p1->neighbors[j].part;
            if(p2->number > i){
                q2 = (bop6 + p2->number * (2 * l + 1) + l);
                calc_harmonics(l, i, p2->number, q1, q2);
            }
        }
        //At this point, bop6 holds the sum of spherical harmonics, a.k.a. ten Wolde's q_{lm}
    }

    //Lechner & Dellago: \bar{q}_{lm}
    for (int i =0; i<N; i++){
        p1 = particles+i;
        int nnneigh1 = p1->nnneigh;
        for (int m=-l; m<=l; m++){
            q1 = (qbar+i * (2*l+1) + m + l);
            q1->re += (*(bop6+i * (2*l+1) + m + l)).re;
            q1->im += (*(bop6+i * (2*l+1) + m + l)).im;
            for(int j = 0; j<nnneigh1; j++){
                p2 = p1->neighbors[j].part;
                q2 = (bop6+p2->number * (2*l+1) + m + l);
                q1->re += q2->re;
                q1->im += q2->im;
            }
            q1->re /= (nnneigh1 + 1.0f);
            q1->im /= (nnneigh1 + 1.0f);
        }
    }

    //Lechner & Dellago: \bar{q}_l
    double prefactor = sqrt(4.0f*M_PI/(2.0f*l+1.0f));
    for (int i=0; i<N; i++){
        temp = sqrt(dotprodsum(qbar+i*(2*l+1), qbar+i*(2*l+1), l));
        (qbar+i*(2*l+1)+l)->re = prefactor*temp;
        
    }
    
    //ten Wolde: \tilde{q}_{lm}(i)
    for(int i = 0; i < N; i++){
        temp = sqrt(dotprodsum(bop6 + i * (2 * l + 1), bop6 + i * (2 * l + 1), l));
        temp = 1.0 / temp;
        for(int m = -l ; m <= l; m++){
            (*(bop6+i * (2 * l + 1) + m + l)).re *= temp;
            (*(bop6+i * (2 * l + 1) + m + l)).im *= temp;
        }
    }
}


/**************************************************
**                  CALC_CONN
** Determines whether a particle is solid-like
**************************************************/
void calc_conn(){
    int z;
    const int l = 6;
    particle *p1;
    particle *p2;
    for(int i = 0; i < N; i++){
        p1 = particles+i;
        z = 0;
        for(int j = 0; j < p1->nnneigh; j++){
            p2 = p1->neighbors[j].part;
            if(dotprodsum(bop6 + i * (2 * l + 1), bop6 + p2->number * (2 * l + 1), l) > 0.7f)
                z++;
        }
        if(z >= 6) p1->type = 1;
        else p1->type = 0;
    }
}


/**************************************************
**                  COMPARE
** Return: 1 if a>b, 0 if a=b, -1 if a<b
**************************************************/
int compare(const void* a, const void* b){
    double d = ((neigh*) a)->dist - ((neigh*) b)->dist;
    return ((0 < d) - (d < 0));
}


/**************************************************
**                  BUILDSANN
** Builds list of nearest neighbors using SANN
**************************************************/
void buildsann(){

    particle *p;
    particle *pj;
    double dx, dy, dz;

    for (int i=0; i<N; i++){
        p = particles+i;
        int numposnb = p->nneigh;
        for (int j=0; j<numposnb; j++){
            pj = p->neighbors[j].part;
            dx = p->x - pj->x;
            dy = p->y - pj->y;
            dz = p->z - pj->z;
            if (dx > 0.5 * xsize) dx -= xsize; else if (dx < -0.5 * xsize) dx += xsize;  //periodic boundaries
            if (dy > 0.5 * ysize) dy -= ysize; else if (dy < -0.5 * ysize) dy += ysize;
            if (dz > 0.5 * zsize) dz -= zsize; else if (dz < -0.5 * zsize) dz += zsize;
            p->neighbors[j].dist = sqrt(dx*dx + dy*dy + dz*dz);
        }

        qsort(p->neighbors, numposnb, sizeof(neigh), compare);

        int m = 3; 
        int done = 0;
        while (!done){
            double rim = 0;
            for (int j = 0; j < m; j++)
                rim += p->neighbors[j].dist / (m - 2);
            if (rim > p->neighbors[m].dist) m++;
            else done = 1;
            if (m > numposnb){ 
                printf("[EDMC-poly] WARNING: NN algorithm did not converge! (m: %d > numposnb: %d)\n", m, numposnb);
                m--;
                done = 1;
            }
        }
        p->nnneigh = m;
    }
}
