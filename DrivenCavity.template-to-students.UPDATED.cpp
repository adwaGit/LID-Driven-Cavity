/************************************************************************ */
/*      This code solves for the viscous flow in a lid-driven cavity      */
/**************************************************************************/

#include <iostream> 
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace std;

/************* Following are fixed parameters for array sizes **************/
#define imax 65     /* Number of points in the x-direction (use odd numbers only) */
#define jmax 65     /* Number of points in the y-direction (use odd numbers only) */
#define neq 3       /* Number of equation to be solved ( = 3: mass, x-mtm, y-mtm) */

/**********************************************/
/****** All Global variables declared here. ***/
/***** These variables SHOULD not be changed **/
/********** by the program once set. **********/
/**********************************************/
/***** The variables declared "const" CAN *****/
/*** not be changed by the program once set ***/
/**********************************************/

/*--------- Numerical Constants --------*/
  const double zero   = 0.0;
  const double tenth  = 0.1;
  const double sixth  = 1.0/6.0;
  const double fifth  = 0.2;
  const double fourth = 0.25;
  const double third  = 1.0/3.0;
  const double half   = 0.5;
  const double one    = 1.0;
  const double two    = 2.0;
  const double three  = 3.0;
  const double four   = 4.0;
  const double six    = 6.0;
  
/*--------- User sets inputs here  --------*/

  const int nmax = 500000;             /* Maximum number of iterations */
  const int iterout = 5000;             /* Number of time steps between solution output */
  const int imms = 1;                   /* Manufactured solution flag: = 1 for manuf. sol., = 0 otherwise */
  const int isgs = 1;                   /* Symmetric Gauss-Seidel  flag: = 1 for SGS, = 0 for point Jacobi */
  const int irstr = 0;                  /* Restart flag: = 1 for restart (file 'restart.in', = 0 for initial run */
  const int ipgorder = 0;               /* Order of pressure gradient: 0 = 2nd, 1 = 3rd (not needed) */
  const int lim = 0;                    /* variable to be used as the limiter sensor (= 0 for pressure) */
  const int residualOut = 10;           /* Number of timesteps between residual output */

  const double cfl  = 0.9;              /* CFL number used to determine time step */
  const double Cx = 0.01;               /* Parameter for 4th order artificial viscosity in x */
  const double Cy = 0.01;               /* Parameter for 4th order artificial viscosity in y */
  const double toler = 1.e-10;          /* Tolerance for iterative residual convergence */
  const double rkappa = 0.1;            /* Time derivative preconditioning constant */
  const double Re = 100.0;              /* Reynolds number = rho*Uinf*L/rmu */
  const double pinf = 0.801333844662;   /* Initial pressure (N/m^2) -> from MMS value at cavity center */
  const double uinf = 1.0;              /* Lid velocity (m/s) */
  const double rho = 1.0;               /* Density (kg/m^3) */
  const double xmin = 0.0;              /* Cavity dimensions...: minimum x location (m) */
  const double xmax = 0.05;             /* maximum x location (m) */
  const double ymin = 0.0;              /* maximum y location (m) */
  const double ymax = 0.05;             /*  maximum y location (m) */
  const double Cx2 = 0.0;               /* Coefficient for 2nd order damping (not required) */
  const double Cy2 = 0.0;               /* Coefficient for 2nd order damping (not required) */
  const double fsmall = 1.e-20;         /* small parameter */

/*-- Derived input quantities (set by function 'set_derived_inputs' called from main)----*/
 
  double rhoinv;    /* Inverse density, 1/rho (m^3/kg) */
  double rlength;   /* Characteristic length (m) [cavity width] */
  double rmu;       /* Viscosity (N*s/m^2) */
  double vel2ref;   /* Reference velocity squared (m^2/s^2) */
  double dx;        /* Delta x (m) */
  double dy;        /* Delta y (m) */
  double rpi;       /* Pi = 3.14159... (defined below) */

/*-- Constants for manufactured solutions ----*/
  const double phi0[neq] = {0.25, 0.3, 0.2};            /* MMS constant */
  const double phix[neq] = {0.5, 0.15, 1.0/6.0};        /* MMS amplitude constant */
  const double phiy[neq] = {0.4, 0.2, 0.25};            /* MMS amplitude constant */
  const double phixy[neq] = {1.0/3.0, 0.25, 0.1};       /* MMS amplitude constant */
  const double apx[neq] = {0.5, 1.0/3.0, 7.0/17.0};     /* MMS frequency constant */
  const double apy[neq] = {0.2, 0.25, 1.0/6.0};         /* MMS frequency constant */
  const double apxy[neq] = {2.0/7.0, 0.4, 1.0/3.0};     /* MMS frequency constant */
  const double fsinx[neq] = {0.0, 1.0, 0.0};            /* MMS constant to determine sine vs. cosine */
  const double fsiny[neq] = {1.0, 0.0, 0.0};            /* MMS constant to determine sine vs. cosine */
  const double fsinxy[neq] = {1.0, 1.0, 0.0};           /* MMS constant to determine sine vs. cosine */
                                                        /* Note: fsin = 1 means the sine function */
                                                        /* Note: fsin = 0 means the cosine function */
                                                        /* Note: arrays here refer to the 3 variables */ 


/*****************************************************************************
*                              Array3 Class
*
*****************************************************************************/

class Array3
{
    private:
        int idim, jdim, kdim;
        double *data;

    public:
    
        Array3(int, int, int);
        ~Array3();

        void copyData(Array3&);
        void swapData(Array3&);     
    
        double& operator() (int, int, int);
        double operator() (int, int, int) const;
};

Array3::Array3 (int i, int j, int k)
{
    idim = i;
    jdim = j;
    kdim = k;
    data = new double[i*j*k];
}

Array3::~Array3 ()
{
    delete [] data;
}

//Copies data from (Array3& A) into the calling Array3 class.   Both Array3's now contain identical data arrays
void Array3::copyData (Array3& A) 
{
    memcpy( data, A.data, idim*jdim*kdim*sizeof(double) );
}


//Swaps pointers to data--thus U.swapData(U2) exchanges data arrays between U and U2
void Array3::swapData (Array3& A)                  
{
    double *temp;

    temp = data;
    data = A.data;
    A.data = temp;
}

inline
double& Array3::operator() (int i, int j, int k)
{
    return data[i*jdim*kdim + j*kdim + k];
    //return data[k*idim*jdim + i*jdim + j];
}

inline      
double Array3::operator() (int i, int j, int k) const
{
    return data[i*jdim*kdim + j*kdim + k];
    //return data[k*idim*jdim + i*jdim + j];
}

/*****************************************************************************
*                              End Array3 Class
*****************************************************************************/



/*****************************************************************************
*                              Array2 Class
*
*****************************************************************************/

class Array2
{
    private:
        int idim, jdim;
        double *data;

    public:
    
        Array2(int, int);
        ~Array2();

        void copyData(Array2&);
        void swapData(Array2&);     
    
        double& operator() (int, int);
        double operator() (int, int) const;
};

Array2::Array2 (int i, int j)
{
    idim = i;
    jdim = j;
    data = new double[i*j];
}

Array2::~Array2 ()
{
    delete [] data;
}

void Array2::copyData (Array2& A)                   //Copies data from (Array2& A) into the calling Array2 class.   
{                                                   //    Both Array2's now contain identical data arrays
    memcpy( data, A.data, idim*jdim*sizeof(double) );
}

void Array2::swapData (Array2& A)                   //Swaps pointers to data--
{                                                   //   thus U.swapData(U2) exchanges data arrays between U and U2
    double *temp;

    temp = data;
    data = A.data;
    A.data = temp;
}

inline
double& Array2::operator() (int i, int j)
{
    return data[i*jdim + j];
}

inline      
double Array2::operator() (int i, int j) const
{
    return data[i*jdim + j];
}

/*****************************************************************************
*                              End Array2 Class
*****************************************************************************/



/*****************Function Pointer Typedefs *********************************/

typedef void (*boundaryConditionPointer)( Array3& );

typedef void (*iterationStepPointer)( boundaryConditionPointer, Array3&, Array3&, Array3&, Array2&, Array2&, Array2& );

/**********************Function Prototypes**********************************/

void set_derived_inputs();
void GS_iteration( boundaryConditionPointer, Array3&, Array3&, Array3&, Array2&, Array2&, Array2& );
void PJ_iteration( boundaryConditionPointer, Array3&, Array3&, Array3&, Array2&, Array2&, Array2& );
void output_file_headers();
void initial( int&, double&, double [neq], Array3&, Array3& );
void bndry( Array3& );
void bndrymms( Array3& );
void write_output( int, Array3&, Array2&, double [neq], double );
double umms( double, double, int ); 
void compute_source_terms( Array3& ); 
double srcmms_mass( double, double );
double srcmms_xmtm( double, double );
double srcmms_ymtm( double, double );
void compute_time_step( Array3&, Array2&, double& );
void Compute_Artificial_Viscosity( Array3&, Array2&, Array2& );
void SGS_forward_sweep( Array3&, Array2&, Array2&, Array2&, Array3& );
void SGS_backward_sweep( Array3&, Array2&, Array2&, Array2&, Array3& );
void point_Jacobi( Array3&, Array3&, Array2&, Array2&, Array2&, Array3& );
void pressure_rescaling( Array3& );
void check_iterative_convergence( int, Array3&, Array3&, Array2&, double [neq], double [neq], int, double, double, double& );
void Discretization_Error_Norms( Array3& );
 

/****************** Inline Function Declarations ***************************/


inline double pow2(double x)                      /* Returns x^2 ... Duplicates pow(x,2)*/
{
    double x2 = x*x;
    return x2;
}

inline double pow3(double x)                      /* Returns x^3 ... Duplicates pow(x,3)*/
{
    double x3 = x*x*x;
    return x3;
}

inline double pow4(double x)                      /* Returns x^4 ... Duplicates pow(x,4)*/
{
    double x4 = x*x*x*x;
    return x4;
}


/******************* End Inline Function Declarations ************************/


/*--- Variables for file handling ---*/
/*--- All files are globally accessible ---*/
  
  FILE *fp1; /* For output of iterative residual history */
  FILE *fp2; /* For output of field data (solution) */
  FILE *fp3; /* For writing the restart file */
  FILE *fp4; /* For reading the restart file */  
  FILE *fp5; /* For output of final DE norms (only for MMS)*/  
//$$$$$$   FILE *fp6; /* For debug: Uncomment for debugging. */  

/***********************************************************************************************************/
/*      NOTE: The Main routine for this C++ code is found at the end                                       */
/***********************************************************************************************************/


/*************************************************************************************************************/
/*                                                                                                           */
/*                                                Functions                                                  */
/*                                                                                                           */
/**********************************************************************************************************  */


void set_derived_inputs()
{
    rhoinv = one/rho;                            /* Inverse density, 1/rho (m^3/kg) */
    rlength = xmax - xmin;                       /* Characteristic length (m) [cavity width] */
    rmu = rho*uinf*rlength/Re;                   /* Viscosity (N*s/m^2) */
    vel2ref = uinf*uinf;                         /* Reference velocity squared (m^2/s^2) */
    dx = (xmax - xmin)/(double)(imax - 1);          /* Delta x (m) */
    dy = (ymax - ymin)/(double)(jmax - 1);          /* Delta y (m) */
    rpi = acos(-one);                            /* Pi = 3.14159... */
    printf("rho,V,L,mu,Re: %f %f %f %f %f\n",rho,uinf,rlength,rmu,Re);
}

/**************************************************************************/

void GS_iteration( boundaryConditionPointer set_boundary_conditions, Array3& u, Array3& uold, Array3& src, Array2& viscx, Array2& viscy, Array2& dt )
{
    /* Copy u to uold (save previous flow values)*/
    uold.copyData(u);

    /* Artificial Viscosity */
    Compute_Artificial_Viscosity(u, viscx, viscy);
              
    /* Symmetric Gauss-Siedel: Forward Sweep */
    SGS_forward_sweep(u, viscx, viscy, dt, src);
          
    /* Set Boundary Conditions for u */
    set_boundary_conditions(u);
           
    /* Artificial Viscosity */
    Compute_Artificial_Viscosity(u, viscx, viscy);
                 
    /* Symmetric Gauss-Siedel: Backward Sweep */
    SGS_backward_sweep(u, viscx, viscy, dt, src);

    /* Set Boundary Conditions for u */
    set_boundary_conditions(u);
}

/**************************************************************************/

void PJ_iteration( boundaryConditionPointer set_boundary_conditions, Array3& u, Array3& uold, Array3& src, Array2& viscx, Array2& viscy, Array2& dt )
{
    /* Swap pointers for u and uold*/
    uold.swapData(u);

    /* Artificial Viscosity */
    Compute_Artificial_Viscosity(uold, viscx, viscy);
              
    /* Point Jacobi: Forward Sweep */
    point_Jacobi(u, uold, viscx, viscy, dt, src);
           
    /* Set Boundary Conditions for u */
    set_boundary_conditions(u);
}

/**************************************************************************/

void output_file_headers()
{
  /*
  Uses global variable(s): imms, fp1, fp2
  */
  
  /* Note: The vector of primitive variables is: */
  /*               u = [p, u, v]^T               */  
  /* Set up output files (history and solution)  */    

    fp1 = fopen("./history.dat","w");
    fprintf(fp1,"TITLE = \"Cavity Iterative Residual History\"\n");
    fprintf(fp1,"variables=\"Iteration\"\"Time(s)\"\"Res1\"\"Res2\"\"Res3\"\n");

    fp2 = fopen("./cavity.dat","w");
    fprintf(fp2,"TITLE = \"Cavity Field Data\"\n");
    if(imms==1)
    {
        fprintf(fp2,"variables=\"x(m)\"\"y(m)\"\"p(N/m^2)\"\"u(m/s)\"\"v(m/s)\"");\
        fprintf(fp2,"\"p-exact\"\"u-exact\"\"v-exact\"\"DE-p\"\"DE-u\"\"DE-v\"\n");      
    }
    else
    {
        if(imms==0)
        {
            fprintf(fp2,"variables=\"x(m)\"\"y(m)\"\"p(N/m^2)\"\"u(m/s)\"\"v(m/s)\"\n");
        }      
        else
        {
            printf("ERROR! imms must equal 0 or 1!!!\n");
            exit (0);
        }       
    }

  /* Header for Screen Output */
  printf("Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum\n"); 
}

/**************************************************************************/

void initial(int& ninit, double& rtime, double resinit[neq], Array3& u, Array3& s)
{
    /* 
    Uses global variable(s): zero, one, irstr, imax, jmax, neq, uinf, pinf 
    To modify: ninit, rtime, resinit, u, s
    */
    int i;                       /* i index (x direction) */
    int j;                       /* j index (y direction) */
    int k;                       /* k index (# of equations) */
  
    double x;       /* Temporary variable for x location */
    double y;       /* Temporary variable for y location */

    /* This subroutine sets inital conditions in the cavity */

    /* Note: The vector of primitive variables is:  */
    /*              u = [p, u, v]^T               */
  
    if(irstr==0)   /* Starting run from scratch*/
    {  
        ninit = 1;          /* set initial iteration to one */
        rtime = zero;       /* set initial time to zero */
        for(k = 0; k<neq; k++)
        {
            resinit[k] = one;
        }
        for(i = 0; i<imax; i++)
        {
            for(j = 0; j<jmax; j++)
            {
                u(i,j,0) = pinf;
                u(i,j,1) = zero;
                u(i,j,2) = zero;

                s(i,j,0) = zero;
                s(i,j,1) = zero;
                s(i,j,2) = zero;
            }
            u(i, jmax-1, 1) = uinf; /* Initialize lid (top) to freestream velocity */
        }
    }  
    else if(irstr==1)  /* Restarting from previous run (file 'restart.in') */
    {
        fp4 = fopen("./restart.in","r"); /* Note: 'restart.in' must exist! */
        if (fp4==NULL)
        {
            printf("Error opening restart file. Stopping.\n");
            exit (0);
        }      
        fscanf(fp4, "%d %lf", ninit, rtime); /* Need to known current iteration # and time value */
        fscanf(fp4, "%lf %lf %lf", &resinit[0], &resinit[1], &resinit[2]); /* Needs initial iterative residuals for scaling */
        for(i=0; i<imax; i++)
        {
            for(j=0; j<jmax; j++)
            {
                fscanf(fp4, "%lf %lf %lf %lf %lf", &x, &y, &u(i,j,0), &u(i,j,1), &u(i,j,2)); 
            }
        }
        ninit += 1;
        printf("Restarting at iteration %d\n", ninit);
        fclose(fp4);
    }   
    else
    {
        printf("ERROR: irstr must equal 0 or 1!\n");
        exit (0);
    }
}


/**************************************************************************/

void bndry( Array3& u )
{
    /* 
    Uses global variable(s): zero, one (not used), two, half, imax, jmax, uinf  
    To modify: u 
    */
    int i;                                          //i index (x direction)
    int j;                                          //j index (y direction)

    /* This applies the cavity boundary conditions */


/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */



}

/**************************************************************************/

void bndrymms( Array3& u )
{
    /* 
    Uses global variable(s): two, imax, jmax, neq, xmax, xmin, ymax, ymin, rlength  
    To modify: u
    */
    int i;                       /* i index (x direction) */
    int j;                       /* j index (y direction) */
    int k;                       /* k index (# of equations) */
  
    double x;       /* Temporary variable for x location */
    double y;       /* Temporary variable for y location */

    /* This applies the cavity boundary conditions for the manufactured solution */

    /* Side Walls */
    for( j = 1; j<jmax-1; j++)
    {
        y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
        i = 0;
        x = xmin;
            
        u(i,j,0)  = umms(x,y,0);
        u(i,j,1)  = umms(x,y,1);
        u(i,j,2)  = umms(x,y,2);

        u(0,j,0) = two*u(1,j,0) - u(2,j,0);    /* 2nd Order BC */

        i=imax-1;
        x = xmax;
            
        u(i,j,0)  = umms(x,y,0);
        u(i,j,1)  = umms(x,y,1);
        u(i,j,2)  = umms(x,y,2);

        u(imax-1,j,0) = two*u(imax-2,j,0) - u(imax-3,j,0);   /* 2nd Order BC */
    }

    /* Top/Bottom Walls */
    for(i=0; i<imax; i++)
    {
        x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
        j = 0;
        y = ymin;

        u(i,j,0)  = umms(x,y,0);
        u(i,j,1)  = umms(x,y,1);
        u(i,j,2)  = umms(x,y,2);

        u(i,0,0) = two*u(i,1,0) - u(i,2,0);   /* 2nd Order BC */

        j = jmax-1;
        y = ymax;
            
        u(i,j,0)  = umms(x,y,0);
        u(i,j,1)  = umms(x,y,1);
        u(i,j,2)  = umms(x,y,2);

        u(i,jmax-1,0) = two*u(i,jmax-2,0) - u(i,jmax-3,0);   /* 2nd Order BC */
    }
}

/**************************************************************************/

void write_output(int n, Array3& u, Array2& dt, double resinit[neq], double rtime)
{
        /* 
    Uses global variable(s): imax, jmax, new, xmax, xmin, ymax, ymin, rlength, imms
    Uses global variable(s): ninit, u, dt, resinit, rtime
    To modify: <none> 
    Writes output and restart files.
    */
   
    int i;                       /* i index (x direction) */
    int j;                       /* j index (y direction) */
    int k;                       /* k index (# of equations) */

    double x;       /* Temporary variable for x location */
    double y;       /* Temporary variable for y location */

    /* Field output */
    fprintf(fp2, "zone T=\"n=%d\"\n",n);
    fprintf(fp2, "I= %d J= %d\n",imax, jmax);
    fprintf(fp2, "DATAPACKING=POINT\n");

    if(imms==1) 
    {
        for(i=0; i<imax; i++)
        {
            for(j=0; j<jmax; j++)
            {
                x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
                y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
                fprintf(fp2,"%e %e %e %e %e %e %e %e %e %e %e\n", x, y, u(i,j,0), u(i,j,1), u(i,j,2), 
                                               umms(x,y,0), umms(x,y,1), umms(x,y,2), 
                                                (u(i,j,0)-umms(x,y,0)), (u(i,j,1)-umms(x,y,1)), (u(i,j,2)-umms(x,y,2)));
            }
        }    
    }
    else if(imms==0)
    {
        for(i=0; i<imax; i++)
        {
            for(j=0; j<jmax; j++)
            {
                x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
                y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
                fprintf(fp2,"%e %e %e %e %e\n", x, y, u(i,j,0), u(i,j,1), u(i,j,2));
            }
        }
    }
    else
    {
        printf("ERROR: imms must equal 0 or 1!\n");
        exit (0);
    }

    /* Restart file: overwrites every 'iterout' iteration */
    fp3 = fopen("./restart.out","w");       
    fprintf(fp3,"%d %e\n", n, rtime);    
    fprintf(fp3,"%e %e %e\n", resinit[0], resinit[1], resinit[2]);
    for(i=0; i<imax; i++)
    {
        for(j=0; j<jmax; j++)
        {
            x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
            y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
            fprintf(fp3,"%e %e %e %e %e\n", x, y, u(i,j,0), u(i,j,1), u(i,j,2));
        }
    }
    fclose(fp3);
}

/**************************************************************************/

double umms(double x, double y, int k)  
{
    /* 
    Uses global variable(s): one, rpi, rlength
    Inputs: x, y, k
    To modify: <none>
    Returns: umms
    */

    double ummstmp; /* Define return value for umms as double precision */

    double termx;      /* Temp variable */
    double termy;      /* Temp variable */
    double termxy;     /* Temp variable */
    double argx;       /* Temp variable */
    double argy;       /* Temp variable */
    double argxy;      /* Temp variable */  
  
    /* This function returns the MMS exact solution */
  
    argx = apx[k]*rpi*x/rlength;
    argy = apy[k]*rpi*y/rlength;
    argxy = apxy[k]*rpi*x*y/rlength/rlength;
    termx = phix[k]*(fsinx[k]*sin(argx)+(one-fsinx[k])*cos(argx));
    termy = phiy[k]*(fsiny[k]*sin(argy)+(one-fsiny[k])*cos(argy));
    termxy = phixy[k]*(fsinxy[k]*sin(argxy)+(one-fsinxy[k])*cos(argxy));
  
    ummstmp = phi0[k] + termx + termy + termxy;
 
    return ummstmp;  
}

/**************************************************************************/

void compute_source_terms( Array3& s )
{
    /* 
    Uses global variable(s): imax, jmax, imms, rlength, xmax, xmin, ymax, ymin
    To modify: s (source terms)
    */

    int i;                       /* i index (x direction) */
    int j;                       /* j index (y direction) */
    
    double x;       /* Temporary variable for x location */
    double y;       /* Temporary variable for y location */

    /* Evaluate Source Terms Once at Beginning (only interior points; will be zero for standard cavity) */
  
    for(i=1; i<imax-1; i++)
    {
        for(j=1; j<jmax-1; j++)
        {
            x = (xmax - xmin)*(double)(i)/(double)(imax - 1);
            y = (ymax - ymin)*(double)(j)/(double)(jmax - 1);
            s(i,j,0) = (double)(imms)*srcmms_mass(x,y);
            s(i,j,1) = (double)(imms)*srcmms_xmtm(x,y);
            s(i,j,2) = (double)(imms)*srcmms_ymtm(x,y);
        }
    }
}

/**************************************************************************/

double srcmms_mass(double x, double y)  
{
    /* 
    Uses global variable(s): rho, rpi, rlength
    Inputs: x, y
    To modify: <none>
    Returns: srcmms_mass
    */

    double srcmasstmp; /* Define return value for srcmms_mass as double precision */

    double dudx;    /* Temp variable: u velocity gradient in x direction */
    double dvdy;  /* Temp variable: v velocity gradient in y direction */
  
    /* This function returns the MMS mass source term */

    dudx = phix[1]*apx[1]*rpi/rlength*cos(apx[1]*rpi*x/rlength) + phixy[1]*apxy[1]*rpi*y/rlength/rlength * cos(apxy[1]*rpi*x*y/rlength/rlength);
  
    dvdy = -phiy[2]*apy[2]*rpi/rlength*sin(apy[2]*rpi*y/rlength) - phixy[2]*apxy[2]*rpi*x/rlength/rlength * sin(apxy[2]*rpi*x*y/rlength/rlength);

    srcmasstmp = rho*dudx + rho*dvdy;

    return srcmasstmp;
} 

/**************************************************************************/

double srcmms_xmtm(double x, double y)  
{
    /* 
    Uses global variable(s): rho, rpi, rmu, rlength
    Inputs: x, y
    To modify: <none>
    Returns: srcmms_xmtm
    */

    double srcxmtmtmp; /* Define return value for srcmms_xmtm as double precision */

    double dudx;    /* Temp variable: u velocity gradient in x direction */
    double dudy;  /* Temp variable: u velocity gradient in y direction */
    double termx;        /* Temp variable */
    double termy;        /* Temp variable */
    double termxy;       /* Temp variable */
    double uvel;         /* Temp variable: u velocity */
    double vvel;         /* Temp variable: v velocity */
    double dpdx;         /* Temp variable: pressure gradient in x direction */
    double d2udx2;       /* Temp variable: 2nd derivative of u velocity in x direction */
    double d2udy2;       /* Temp variable: 2nd derivative of u velocity in y direction */

    /*This function returns the MMS x-momentum source term */

    termx = phix[1]*sin(apx[1]*rpi*x/rlength);
    termy = phiy[1]*cos(apy[1]*rpi*y/rlength);
    termxy = phixy[1]*sin(apxy[1]*rpi*x*y/rlength/rlength);
    uvel = phi0[1] + termx + termy + termxy;
  
    termx = phix[2]*cos(apx[2]*rpi*x/rlength);
    termy = phiy[2]*cos(apy[2]*rpi*y/rlength);
    termxy = phixy[2]*cos(apxy[2]*rpi*x*y/rlength/rlength);
    vvel = phi0[2] + termx + termy + termxy;
  
    dudx = phix[1]*apx[1]*rpi/rlength*cos(apx[1]*rpi*x/rlength) + phixy[1]*apxy[1]*rpi*y/rlength/rlength * cos(apxy[1]*rpi*x*y/rlength/rlength);
  
    dudy = -phiy[1]*apy[1]*rpi/rlength*sin(apy[1]*rpi*y/rlength) + phixy[1]*apxy[1]*rpi*x/rlength/rlength * cos(apxy[1]*rpi*x*y/rlength/rlength);
  
    dpdx = -phix[0]*apx[0]*rpi/rlength*sin(apx[0]*rpi*x/rlength) + phixy[0]*apxy[0]*rpi*y/rlength/rlength * cos(apxy[0]*rpi*x*y/rlength/rlength);

    d2udx2 = -phix[1]*pow((apx[1]*rpi/rlength),2) * sin(apx[1]*rpi*x/rlength) - phixy[1]*pow((apxy[1]*rpi*y/rlength/rlength),2) * sin(apxy[1]*rpi*x*y/rlength/rlength); 
 
    d2udy2 = -phiy[1]*pow((apy[1]*rpi/rlength),2) * cos(apy[1]*rpi*y/rlength) - phixy[1]*pow((apxy[1]*rpi*x/rlength/rlength),2) * sin(apxy[1]*rpi*x*y/rlength/rlength);
  
    srcxmtmtmp = rho*uvel*dudx + rho*vvel*dudy + dpdx - rmu*( d2udx2 + d2udy2 );

    return srcxmtmtmp;
} 

/**************************************************************************/

double srcmms_ymtm(double x, double y)  
{
    /* 
    Uses global variable(s): rho, rpi, rmu, rlength
    Inputs: x, y
    To modify: <none>
    Returns: srcmms_ymtm
    */

    double srcymtmtmp; /* Define return value for srcmms_ymtm as double precision */

    double dvdx;         /* Temp variable: v velocity gradient in x direction */
    double dvdy;         /* Temp variable: v velocity gradient in y direction */
    double termx;        /* Temp variable */
    double termy;        /* Temp variable */
    double termxy;       /* Temp variable */
    double uvel;         /* Temp variable: u velocity */
    double vvel;         /* Temp variable: v velocity */
    double dpdy;         /* Temp variable: pressure gradient in y direction */
    double d2vdx2;       /* Temp variable: 2nd derivative of v velocity in x direction */
    double d2vdy2;       /* Temp variable: 2nd derivative of v velocity in y direction */

    /* This function returns the MMS y-momentum source term */

    termx = phix[1]*sin(apx[1]*rpi*x/rlength);
    termy = phiy[1]*cos(apy[1]*rpi*y/rlength);
    termxy = phixy[1]*sin(apxy[1]*rpi*x*y/rlength/rlength);
    uvel = phi0[1] + termx + termy + termxy;
  
    termx = phix[2]*cos(apx[2]*rpi*x/rlength);
    termy = phiy[2]*cos(apy[2]*rpi*y/rlength);
    termxy = phixy[2]*cos(apxy[2]*rpi*x*y/rlength/rlength);
    vvel = phi0[2] + termx + termy + termxy;
  
    dvdx = -phix[2]*apx[2]*rpi/rlength*sin(apx[2]*rpi*x/rlength) - phixy[2]*apxy[2]*rpi*y/rlength/rlength * sin(apxy[2]*rpi*x*y/rlength/rlength);
  
    dvdy = -phiy[2]*apy[2]*rpi/rlength*sin(apy[2]*rpi*y/rlength) - phixy[2]*apxy[2]*rpi*x/rlength/rlength * sin(apxy[2]*rpi*x*y/rlength/rlength);
  
    dpdy = phiy[0]*apy[0]*rpi/rlength*cos(apy[0]*rpi*y/rlength) + phixy[0]*apxy[0]*rpi*x/rlength/rlength * cos(apxy[0]*rpi*x*y/rlength/rlength);
  
    d2vdx2 = -phix[2]*pow((apx[2]*rpi/rlength),2) * cos(apx[2]*rpi*x/rlength) - phixy[2]*pow((apxy[2]*rpi*y/rlength/rlength),2) * cos(apxy[2]*rpi*x*y/rlength/rlength);
  
    d2vdy2 = -phiy[2]*pow((apy[2]*rpi/rlength),2) * cos(apy[2]*rpi*y/rlength)   - phixy[2]*pow((apxy[2]*rpi*x/rlength/rlength),2) * cos(apxy[2]*rpi*x*y/rlength/rlength);
  
    srcymtmtmp = rho*uvel*dvdx + rho*vvel*dvdy + dpdy - rmu*( d2vdx2 + d2vdy2 );
  
    return srcymtmtmp;  
}

/**************************************************************************/

void compute_time_step( Array3& u, Array2& dt, double& dtmin )
{
    /* 
    Uses global variable(s): one (not used), two, four, half, fourth
    Uses global variable(s): vel2ref, rmu, rho, dx, dy, cfl, rkappa, imax, jmax
    Uses: u
    To Modify: dt, dtmin
    */
    int i;                      //i index (x direction)
    int j;                      //j index (y direction)
  
    double dtvisc;          //Viscous time step stability criteria (constant over domain)
    double uvel2;           //Local velocity squared
    double beta2;           //Beta squared parameter for time derivative preconditioning
    double lambda_x;        //Max absolute value eigenvalue in (x,t)
    double lambda_y;        //Max absolute value eigenvalue in (y,t)
    double lambda_max;      //Max absolute value eigenvalue (used in convective time step computation)
    double dtconv;          //Local convective time step restriction

/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */



}  

/**************************************************************************/

void Compute_Artificial_Viscosity( Array3& u, Array2& viscx, Array2& viscy )
{
    /* 
    Uses global variable(s): zero (not used), one (not used), two, four, six, half, fourth (not used)
    Uses global variable(s): imax, jmax, lim (not used), rho, dx, dy, Cx, Cy, Cx2 (not used), Cy2 (not used)
    , fsmall (not used), vel2ref, rkappa
    Uses: u
    To Modify: artviscx, artviscy
    */
    int i;                  //i index (x direction)
    int j;                  //j index (y direction)

    double uvel2;       //Local velocity squared
    double beta2;       //Beta squared parameter for time derivative preconditioning
    double lambda_x;    //Max absolute value e-value in (x,t)
    double lambda_y;    //Max absolute value e-value in (y,t)
    double d4pdx4;      //4th derivative of pressure w.r.t. x
    double d4pdy4;      //4th derivative of pressure w.r.t. y

/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */



}

/**************************************************************************/

void SGS_forward_sweep( Array3& u, Array2& viscx, Array2& viscy, Array2& dt, Array3& s )
{
    /* 
    Uses global variable(s): two, three (not used), six (not used), half
    Uses global variable(s): imax, jmax, ipgorder (not used), rho, rhoinv, dx, dy, rkappa,
                        xmax (not used), xmin (not used), ymax (not used), ymin (not used), 
                        rmu, vel2ref
    Uses: artviscx, artviscy, dt, s
    To Modify: u
    */
 
    double dpdx;        //First derivative of pressure w.r.t. x
    double dudx;        //First derivative of x velocity w.r.t. x
    double dvdx;        //First derivative of y velocity w.r.t. x
    double dpdy;        //First derivative of pressure w.r.t. y
    double dudy;        //First derivative of x velocity w.r.t. y
    double dvdy;        //First derivative of y velocity w.r.t. y
    double d2udx2;      //Second derivative of x velocity w.r.t. x
    double d2vdx2;      //Second derivative of y velocity w.r.t. x
    double d2udy2;      //Second derivative of x velocity w.r.t. y
    double d2vdy2;      //Second derivative of y velocity w.r.t. y
    double uvel2;       //velocity squared at node
    double beta2;       //Beta squared parameter for time derivative preconditioning

    /* Symmetric Gauss-Siedel: Forward Sweep */

/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */




}

/**************************************************************************/

void SGS_backward_sweep( Array3& u, Array2& viscx, Array2& viscy, Array2& dt, Array3& s )
{
    /* 
    Uses global variable(s): two, three (not used), six (not used), half
    Uses global variable(s): imax, jmax, ipgorder (not used), rho, rhoinv, dx, dy, rkappa,
                        xmax (not used), xmin (not used), ymax (not used), ymin (not used), 
                        rmu, vel2ref
    Uses: artviscx, artviscy, dt, s
    To Modify: u
    */
 
    double dpdx;        //First derivative of pressure w.r.t. x
    double dudx;        //First derivative of x velocity w.r.t. x
    double dvdx;        //First derivative of y velocity w.r.t. x
    double dpdy;        //First derivative of pressure w.r.t. y
    double dudy;        //First derivative of x velocity w.r.t. y
    double dvdy;        //First derivative of y velocity w.r.t. y
    double d2udx2;      //Second derivative of x velocity w.r.t. x
    double d2vdx2;      //Second derivative of y velocity w.r.t. x
    double d2udy2;      //Second derivative of x velocity w.r.t. y
    double d2vdy2;      //Second derivative of y velocity w.r.t. y
    double uvel2;       //Velocity squared at node
    double beta2;       //Beta squared parameter for time derivative preconditioning

    /* Symmetric Gauss-Siedel: Backward Sweep  */


/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */




}

/**************************************************************************/

void point_Jacobi( Array3& u, Array3& uold, Array2& viscx, Array2& viscy, Array2& dt, Array3& s )
{
    /* 
    Uses global variable(s): two, three (not used), six (not used), half
    Uses global variable(s): imax, jmax, ipgorder (not used), rho, rhoinv, dx, dy, rkappa,
                        xmax (not used), xmin (not used), ymax (not used), ymin (not used), 
                        rmu, vel2ref
    Uses: uold, artviscx, artviscy, dt, s
    To Modify: u
    */
 
    double dpdx;        //First derivative of pressure w.r.t. x
    double dudx;        //First derivative of x velocity w.r.t. x
    double dvdx;        //First derivative of y velocity w.r.t. x
    double dpdy;        //First derivative of pressure w.r.t. y
    double dudy;        //First derivative of x velocity w.r.t. y
    double dvdy;        //First derivative of y velocity w.r.t. y
    double d2udx2;      //Second derivative of x velocity w.r.t. x
    double d2vdx2;      //Second derivative of y velocity w.r.t. x
    double d2udy2;      //Second derivative of x velocity w.r.t. y
    double d2vdy2;      //Second derivative of y velocity w.r.t. y
    double uvel2;       //Velocity squared at node
    double beta2;       //Beta squared parameter for time derivative preconditioning

    /* Point Jacobi method */


/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */




}

/**************************************************************************/

void pressure_rescaling( Array3& u )
{
    /* 
    Uses global variable(s): imax, jmax, imms, xmax, xmin, ymax, ymin, rlength, pinf
    To Modify: u
    */

    int iref;                     /* i index location of pressure rescaling point */
    int jref;                     /* j index location of pressure rescaling point */

    double x;               /* Temporary variable for x location */
    double y;               /* Temporary variable for y location */  
    double deltap;          /* delta_pressure for rescaling all values */

    iref = (imax-1)/2;     /* Set reference pressure to center of cavity */
    jref = (jmax-1)/2;
    if(imms==1)
    {
        x = (xmax - xmin)*(double)(iref)/(double)(imax - 1);
        y = (ymax - ymin)*(double)(jref)/(double)(jmax - 1);
        deltap = u(iref,jref,0) - umms(x,y,0); /* Constant in MMS */
    }
    else
    {
        deltap = u(iref,jref,0) - pinf; /* Reference pressure */
    }

    for(int i=0; i<imax; i++)
    {
        for(int j=0; j<jmax; j++)
        {
            u(i,j,0) -= deltap;
        }
    }      
}  

/**************************************************************************/

void check_iterative_convergence(int n, Array3& u, Array3& uold, Array2& dt, double res[neq], double resinit[neq], int ninit, double rtime, double dtmin, double& conv)
{
  /* 
  Uses global variable(s): zero
  Uses global variable(s): imax, jmax, neq, fsmall (not used)
  Uses: n, u, uold, dt, res, resinit, ninit, rtime, dtmin
  To modify: conv
  */

  int i;                       /* i index (x direction) */
  int j;                       /* j index (y direction) */
  int k;                       /* k index (# of equations) */

  /* Compute iterative residuals to monitor iterative convergence */

    res[0] = zero;              //Reset to zero (as they are sums)
    res[1] = zero;
    res[2] = zero;

/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */




    /* Write iterative residuals every "residualOut" iterations */
    if( ((n%residualOut)==0)||(n==ninit) )
    {
        fprintf(fp1, "%d %e %e %e %e\n",n, rtime, res[0], res[1], res[2] );
        printf("%d   %e   %e   %e   %e   %e\n",n, rtime, dtmin, res[0], res[1], res[2] );    

        /* Write header for iterative residuals every 20 residual printouts */
        if( ((n%(residualOut*20))==0)||(n==ninit) )
        {
            printf("Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum\n"); 
        }    
    }
     
}

/**************************************************************************/

void Discretization_Error_Norms( Array3& u ) 
{
    /* 
    Uses global variable(s): zero
    Uses global variable(s): imax, jmax, neq, imms, xmax, xmin, ymax, ymin, rlength (not used)
    Uses: u
    To modify: rL1norm, rL2norm, rLinfnorm 
    */

    double x;                   //Temporary variable for x location
    double y;                   //Temporary variable for y location
    double DE;                      //Discretization error (absolute value)

    double rL1norm[neq];
    double rL2norm[neq]; 
    double rLinfnorm[neq];

    /* Only compute discretization error norms for manufactured solution */
    if(imms==1)
    {

/* !************************************************************** */
/* !************ADD CODING HERE FOR INTRO CFD STUDENTS************ */
/* !************************************************************** */




    }
}

/********************************************************************************************************************/
/*                                                                                                                  */
/*                                                End Functions                                                     */
/*                                                                                                                  */
/********************************************************************************************************************/







/********************************************************************************************************************/
/*                                                                                                                  */
/*                                                Main Function                                                     */
/*                                                                                                                  */
/********************************************************************************************************************/
int main()
{
    //Data class declarations: hold all the data needed across the entire grid
    Array3 u     (imax, jmax, neq);     //u and uold store the current and previous primitive variable solution on the entire grid
    Array3 uold  (imax, jmax, neq);

    Array3 src   (imax, jmax, neq);     //src stores the source terms over the entire grid (used for MMS)

    Array2 viscx (imax, jmax);          //Artificial viscosity, x and y directions
    Array2 viscy (imax, jmax);

    Array2 dt    (imax, jmax);          //Local timestep array



    /* Minimum of iterative residual norms from three equations */
    double conv;
    double resTest;
    int n = 0;  //Iteration number

                                                      
    /*--------- Solution variables declaration ----------------------*/
      
     int ninit = 0;                 /* Initial iteration number (used for restart file) */

     double res[neq];               /* Iterative residual for each equation */
     double resinit[neq];           /* Initial iterative residual for each equation (from iteration 1) */
     double rtime;                  /* Variable to estimate simulation time */
     double dtmin = 1.0e99;         /* Minimum time step for a given iteration (initialized large) */

     double x;                      /* Temporary variable for x location */
     double y;                      /* Temporary variable for y location */


    /*-------Set Function Pointers-----------------------------------*/
    
    iterationStepPointer     iterationStep;
    boundaryConditionPointer set_boundary_conditions;

    if(isgs==1)                 /* ==Symmetric Gauss Seidel== */
    {
        iterationStep = &GS_iteration;
    }
    else if(isgs==0)             /* ==Point Jacobi== */
    {
        iterationStep = &PJ_iteration;
    }
    else
    {
        printf("ERROR: isgs must equal 0 or 1!\n");
        exit (0);  
    }
      
    if(imms==0) 
    {
            set_boundary_conditions = &bndry;
    }
    else if(imms==1)
        {
            set_boundary_conditions = &bndrymms;
        }
        else
        {
            printf("ERROR: imms must equal 0 or 1!\n");
            exit (0);
        }

    /*-------End Set Function Pointers-------------------------------*/

    /* Debug output: Uncomment and modify if debugging */
    //$$$$$$ fp6 = fopen("./Debug.dat","w");
    //$$$$$$ fprintf(fp6,"TITLE = \"Debug Data Data\"\n");
    //$$$$$$ fprintf(fp6,"variables=\"x(m)\"\"y(m)\"\"visc-x\"\"visc-y\"\n");
    //$$$$$$ fprintf(fp6, "zone T=\"n=%d\"\n",n);
    //$$$$$$ fprintf(fp6, "I= %d J= %d\n",imax, jmax);
    //$$$$$$ fprintf(fp6, "DATAPACKING=POINT\n");

    /* Set derived input quantities */
    set_derived_inputs();

    /* Set up headers for output files */
    output_file_headers();

    /* Set Initial Profile for u vector */
    initial( ninit, rtime, resinit, u, src );   

    /* Set Boundary Conditions for u */
    set_boundary_conditions( u );

    /* Write out inital conditions to solution file */
    write_output(ninit, u, dt, resinit, rtime);
     
    /* Evaluate Source Terms Once at Beginning */
    /*(only interior points; will be zero for standard cavity) */
    compute_source_terms( src );

    /*========== Main Loop ==========*/
    for (n = ninit; n<= nmax; n++)
    {
        /* Calculate time step */  
        compute_time_step( u, dt, dtmin );
           
        /* Perform main iteration step (point jacobi or gauss seidel)*/    
        iterationStep( set_boundary_conditions, u, uold, src, viscx, viscy, dt ); 

        /* Pressure Rescaling (based on center point) */
        pressure_rescaling( u );

        /* Update the time */
        rtime += dtmin;

        /* Check iterative convergence using L2 norms of iterative residuals */
        check_iterative_convergence(n, u, uold, dt, res, resinit, ninit, rtime, dtmin, conv);

        if(conv<toler) 
        {
            fprintf(fp1, "%d %e %e %e %e\n",n, rtime, res[0], res[1], res[2]);
                goto converged;
        }
            
        /* Output solution and restart file every 'iterout' steps */
        if( ((n%iterout)==0) ) 
        {
                write_output(n, u, dt, resinit, rtime);
        }
        
    }  /* ========== End Main Loop ========== */

    printf("\nSolver stopped in %d iterations because the specified maximum number of timesteps was exceeded.\n", nmax);
        
    goto notconverged;
        
converged:  /* go here once solution is converged */

    printf("\nSolver stopped in %d iterations because the convergence criteria was met OR because the solution diverged.\n", n);
    printf("   Solution divergence is indicated by inf or NaN residuals.\n", n);
    
notconverged:

    /* Calculate and Write Out Discretization Error Norms (will do this for MMS only) */
    Discretization_Error_Norms( u );

    /* Output solution and restart file */
    write_output(n, u, dt, resinit, rtime);

    /* Close open files */
    fclose(fp1);
    fclose(fp2);
    //$$$$$$   fclose(fp6); /* Uncomment for debug output */

    return 0;
}

