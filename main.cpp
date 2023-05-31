#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

//declare global variables for use in each function
const int Nx=17;          //number of rows
const int Ny=34;          //number of columns
const double g = 9.81;    //gravity constant near earth's surface
const double F0=0.0001;   //Coriolis parameter constant for a latitude of 45 degrees where the calculations occur
double SF[Nx][Ny];        //array for the stream field or psi at time zer0: g/f0*height field
double HT[Nx][Ny];        //height field array, read in from file msc11.dat
double Z[Nx][Ny];         //zeta function array
double psi_py[Nx][Ny];    //U function array, holding -psi_y for each point in the matrix
double psi_px[Nx][Ny];    //V function array, holding psi_x for each point in the matrix
double EDDY[Nx][Ny];      //Eddy viscosity dissipation array holding the values for kappa* det^2(zeta)
double DX;                //distance between matrix points in km
double KAPPA=0.0001;      //eddy diffusivity
double BETA=(sqrt(2)*.00007292)/6370000;    //beta parameter which is equal to omega*sin(phi)/a
                                                   //where omega is earth's rotational speed, phi is the latitude and
                                                   //a is earth's radius
double zeta_pt[Nx][Ny];   //array holding the RHS of the barotropic vorticity eq


void HF();          //height field function
void PSI();         //initial stream field function
void ZETA();        //det^2(psi) function
void U();           //negative partial y derivative of psi
void V();           //partial x derivative of psi
void EVD();         //eddy viscosity dissipation fucntion
void F();           //RHS of barotropic vorticity equation or the value of psi_t
void Z_K(double DT);       //zeta(k+1)
void print_Z();            //prints the matrix zeta


int main ()
{

    //run each function once to calculate everything at time zero
    HF();
    PSI();
    ZETA();
    U();
    V();
    EVD();
    F();

    cout<<scientific;       //puts the output in scientific notation
    //int a=32;             //when using ofstream, a can be used to add a different symbol at the end of each file name

    //char d[]="matrix";

    for(int i=0; i<=48; i++)        //loop runs a forecast using .5 as DT and for a 24 hour period
    {
        cout<<"\ntime: "<<i*.5<<'\n';
        print_Z();
        Z_K(.5);
        EVD();
        F();
       // d[6]=static_cast<char>(a++);
    }


    return 0;
}

/*
 * reads in the height field from file msc11.dat, as well as the row, col, and delta x (which is the same as delta y)
 * height field is in the form of a 17x34 matrix
 * distance between points in physical terms is DX and is 3.125 km
 */
void HF()
{

    int row,col;
    ifstream is;
    is.open("msc11.dat");
    is>>col;
    is>>row;
    is>>DX;

    for(int i=0;i<Nx;i++)     //first line contains parameter sizes and formats, so data starts on second line
    {
        is.ignore();
        for(int j=0; j<Ny; j++)
        {
            is>>HT[i][j];
        }
    }

}

/*
 * Function computes the stream field at time zero using the height field HT and the gravity constant 9.81 and Coriolis
 * parameter F0;
 */
void PSI()
{
    for(int i=0; i<Nx; i++)
    {
        for(int j=0; j<Ny; j++)
        {
            SF[i][j]=(HT[i][j]) * (g/F0);
        }
    }
}

/*
 * Uses det squared psi to approximate a value for the relative velocity
 * Forward differencing is used when ever on the first row or column, backward differencing is used when on the last
 * row and column and centered differencing is used for all points inside the mesh
 *
 * Finds the numerical approximation for the sum of a second order derivative of psi wrt x and second order derivative
 * of psi wrt y:
 *
 * partial^2 psi/partial x^2  +  partial^2 psi/partial y^2
 */
void ZETA()
{

    for(int i=0; i<Nx; i++)
    {
        for(int j=0; j<Ny; j++){
            if(i==0 && j<Ny-3)
                Z[i][j]=(4*(SF[i][j]+SF[i+2][j]+SF[i][j+2])-5*(SF[i+1][j]+SF[i][j+1])-(SF[i+3][j]+SF[i][j+3]))/(pow(DX,3));  //forward diff
            else if(i==0 && j>=Ny-3)
                Z[i][j]=(4*(SF[i][j]+SF[i-2][j]+SF[i][j-2])-5*(SF[i-1][j]+SF[i][j-1])-(SF[i-3][j]+SF[i][j-3]))/(pow(DX,3));  //backward diff
            else if(j==0 && i<Nx-3)
                Z[i][j]=(4*(SF[i][j]+SF[i+2][j]+SF[i][j+2])-5*(SF[i+1][j]+SF[i][j+1])-(SF[i+3][j]+SF[i][j+3]))/(pow(DX,3));  //forward diff
            else if(j==0 && i>=Nx-3)
                Z[i][j]=(4*(SF[i][j]+SF[i-2][j]+SF[i][j-2])-5*(SF[i-1][j]+SF[i][j-1])-(SF[i-3][j]+SF[i][j-3]))/(pow(DX,3));  //backward diff
            else if(i==Nx-1 && j<Ny-3)
                Z[i][j]=(4*(SF[i][j]+SF[i+2][j]+SF[i][j+2])-5*(SF[i+1][j]+SF[i][j+1])-(SF[i+3][j]+SF[i][j+3]))/(pow(DX,3));  //forward diff
            else if(i==Nx-1 && j>=Ny-3)
                Z[i][j]=(4*(SF[i][j]+SF[i-2][j]+SF[i][j-2])-5*(SF[i-1][j]+SF[i][j-1])-(SF[i-3][j]+SF[i][j-3]))/(pow(DX,3));  //backward
            else if(j==Ny-1 && i<Nx-3)
                Z[i][j]=(4*(SF[i][j]+SF[i+2][j]+SF[i][j+2])-5*(SF[i+1][j]+SF[i][j+1])-(SF[i+3][j]+SF[i][j+3]))/(pow(DX,3));  //forward diff
            else if(j==Ny-1 && i>=Nx-3)
                Z[i][j]=(4*(SF[i][j]+SF[i-2][j]+SF[i][j-2])-5*(SF[i-1][j]+SF[i][j-1])-(SF[i-3][j]+SF[i][j-3]))/(pow(DX,3));  //backward diff
            else
                Z[i][j]=(SF[i+1][j]+SF[i-1][j]+SF[i][j+1]+SF[i][j-1]-(4*SF[i][j]))/(pow(DX,2));                     //centered differencing
        }
    }
}

/*
 * Uses the same forward, backward and centered differencing as the ZETA() function above
 *
 * Approximates the value of -(partial psi / partial y)
 */

void U()
{
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if(i==0 && j<Ny-2)
                psi_py[i][j]=-((-3*SF[i][j])+(4*SF[i][j+1])-SF[i][j+2])/(2*DX);  //forward diff
            else if(i==0 && j>=Ny-2)
                psi_py[i][j]=-((3*SF[i][j])-(4*SF[i][j-1])+SF[i][j-2])/(2*DX);   //backward diff
            else if(j==0 && i<Nx-2)
                psi_py[i][j]=-((-3*SF[i][j])+(4*SF[i][j+1])-SF[i][j+2])/(2*DX);  //forward diff
            else if(j==0 && i>=Ny-2)
                psi_py[i][j]=-((3*SF[i][j])-(4*SF[i][j-1])+SF[i][j-2])/(2*DX);   //backward diff
            else if(i==Nx-1 && j<Ny-2)
                psi_py[i][j]=-((-3*SF[i][j])+(4*SF[i][j+1])-SF[i][j+2])/(2*DX);  //forward diff
            else if(i==Nx-1 && j>=Ny-2)
                psi_py[i][j]=-((3*SF[i][j])-(4*SF[i][j-1])+SF[i][j-2])/(2*DX);   //backward diff
            else if(j==Ny-1 && i<Nx-2)
                psi_py[i][j]=-((-3*SF[i][j])+(4*SF[i][j+1])-SF[i][j+2])/(2*DX);  //forward diff
            else if(j==Ny-1 && i>=Nx-2)
                psi_py[i][j]=-((3*SF[i][j])-(4*SF[i][j-1])+SF[i][j-2])/(2*DX);   //forward diff
            else
                psi_py[i][j]=-(SF[i][j+1]-SF[i][j-1])/(2*DX);                    //centered diff
        }
    }
}

/*
 * Same forward, backward, centered differencing
 *
 * Approximates the value of (partial psi / partial x)
 */
void V()
{
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if(i==0 && j<Ny-2)
                psi_px[i][j]=((-3*SF[i][j])+(4*SF[i+1][j])-SF[i+2][j])/(2*DX);
            else if(i==0 && j>=Ny-2)
                psi_px[i][j]=((3*SF[i][j])-(4*SF[i-1][j])+SF[i-2][j])/(2*DX);
            else if(j==0 && i<Nx-2)
                psi_px[i][j]=((-3*SF[i][j])+(4*SF[i+1][j])-SF[i+2][j])/(2*DX);
            else if(j==0 && i>=Ny-2)
                psi_px[i][j]=((3*SF[i][j])-(4*SF[i-1][j])+SF[i-2][j])/(2*DX);
            else if(i==Nx-1 && j<Ny-2)
                psi_px[i][j]=((-3*SF[i][j])+(4*SF[i+1][j])-SF[i+2][j])/(2*DX);
            else if(i==Nx-1 && j>=Ny-2)
                psi_px[i][j]=((3*SF[i][j])-(4*SF[i-1][j])+SF[i-2][j])/(2*DX);
            else if(j==Ny-1 && i<Nx-2)
                psi_px[i][j]=((-3*SF[i][j])+(4*SF[i+1][j])-SF[i+2][j])/(2*DX);
            else if(j==Ny-1 && i>=Nx-2)
                psi_px[i][j]=((3*SF[i][j])-(4*SF[i-1][j])+SF[i-2][j])/(2*DX);
            else
                psi_px[i][j]=(SF[i+1][j]-SF[i-1][j])/(2*DX);
        }
    }

}

/*
 * Eddy Viscosity Dissipation:
 * The eddy viscosity dissipation with eddy diffusivity kappa is included to prevent excessive generation of small-scale
 * noise caused by the non-linear computational instability.
 *
 * EVD is defined as kappa * det^2(zeta)
 *
 * calculated using forward, backward and centered differencing of the sum of zeta_xx and zeta_yy (where subscript xx
 * and yy notate double partial derivatives).
 */
void EVD()
{
    for(int i=0; i<Nx; i++)
    {
        for(int j=0; j<Ny; j++){
            if(i==0 && j<Ny-3)
                EDDY[i][j]=(4*(Z[i][j]+Z[i+2][j]+Z[i][j+2])-5*(Z[i+1][j]+Z[i][j+1])-(Z[i+3][j]+Z[i][j+3]))/(pow(DX,3));  //forward diff
            else if(i==0 && j>=Ny-3)
                EDDY[i][j]=(4*(Z[i][j]+Z[i-2][j]+Z[i][j-2])-5*(Z[i-1][j]+Z[i][j-1])-(Z[i-3][j]+Z[i][j-3]))/(pow(DX,3));  //backward diff
            else if(j==0 && i<Nx-3)
                EDDY[i][j]=(4*(Z[i][j]+Z[i+2][j]+Z[i][j+2])-5*(Z[i+1][j]+Z[i][j+1])-(Z[i+3][j]+Z[i][j+3]))/(pow(DX,3));  //forward diff
            else if(j==0 && i>=Nx-3)
                EDDY[i][j]=(4*(Z[i][j]+Z[i-2][j]+Z[i][j-2])-5*(Z[i-1][j]+Z[i][j-1])-(Z[i-3][j]+Z[i][j-3]))/(pow(DX,3));  //backward diff
            else if(i==Nx-1 && j<Ny-3)
                EDDY[i][j]=(4*(Z[i][j]+Z[i+2][j]+Z[i][j+2])-5*(Z[i+1][j]+Z[i][j+1])-(Z[i+3][j]+Z[i][j+3]))/(pow(DX,3));  //forward diff
            else if(i==Nx-1 && j>=Ny-3)
                EDDY[i][j]=(4*(Z[i][j]+Z[i-2][j]+Z[i][j-2])-5*(Z[i-1][j]+Z[i][j-1])-(Z[i-3][j]+Z[i][j-3]))/(pow(DX,3));  //backward diff
            else if(j==Ny-1 && i<Nx-3)
                EDDY[i][j]=(4*(Z[i][j]+Z[i+2][j]+Z[i][j+2])-5*(Z[i+1][j]+Z[i][j+1])-(Z[i+3][j]+Z[i][j+3]))/(pow(DX,3));  //forward diff
            else if(j==Ny-1 && i>=Nx-3)
                EDDY[i][j]=(4*(Z[i][j]+Z[i-2][j]+Z[i][j-2])-5*(Z[i-1][j]+Z[i][j-1])-(Z[i-3][j]+Z[i][j-3]))/(pow(DX,3));  //backward diff
            else
                EDDY[i][j]=(Z[i+1][j]+Z[i-1][j]+Z[i][j+1]+Z[i][j-1]-(4*Z[i][j]))/(pow(DX,2));                     //centered differencing
        }
    }

}

/*
 * F() is a function that calculates the entirety of the RHS of the barotropic vorticity equation:
 *      -((u*zeta)_x + (v*zeta)_y + (beta * v)) + (kappa * det^2(zeta))
 *
 * Uses forward, backward and centered differencing to find the first partial derivative in x of (u*zeta) and the
 * first partial derivative in y of (v*zeta)
 *
 * Then adds the (i,j) value of beta * v and the EVD to compute zeta_t
 */

void F()
{

    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if(i==0 && j<Ny-1)             //forward
                zeta_pt[i][j]=-(((-3*psi_py[i][j]*Z[i][j])+(4*psi_py[i+1][j]*Z[i+1][j])-(psi_py[i+2][j]*Z[i+2][j]))+((-3*psi_px[i][j]*Z[i][j])+(4*psi_px[i][j+1]*Z[i][j+1])-(psi_px[i][j+2]*Z[i][j+2])))/(2*DX)-(BETA*psi_px[i][j])+(KAPPA*EDDY[i][j]);
            else if(i==0 && j==Ny-1)       //backward
                zeta_pt[i][j]=-(((3*psi_py[i][j]*Z[i][j])-(4*psi_py[i-1][j]*Z[i-1][j])+(psi_py[i-2][j]*Z[i-2][j]))+((3*psi_px[i][j]*Z[i][j])-(4*psi_px[i][j-1]*Z[i][j-1])+(psi_px[i][j-2]*Z[i][j-2])))/(2*DX)-(BETA*psi_px[i][j])+(KAPPA*EDDY[i][j]);
            else if(j==0 && i<Nx-1)        //forward
                zeta_pt[i][j]=-(((-3*psi_py[i][j]*Z[i][j])+(4*psi_py[i+1][j]*Z[i+1][j])-(psi_py[i+2][j]*Z[i+2][j]))+((-3*psi_px[i][j]*Z[i][j])+(4*psi_px[i][j+1]*Z[i][j+1])-(psi_px[i][j+2]*Z[i][j+2])))/(2*DX)-(BETA*psi_px[i][j])+(KAPPA*EDDY[i][j]);
            else if(j==0 && i==Nx-1)       //backward
                zeta_pt[i][j]=-(((3*psi_py[i][j]*Z[i][j])-(4*psi_py[i-1][j]*Z[i-1][j])+(psi_py[i-2][j]*Z[i-2][j]))+((3*psi_px[i][j]*Z[i][j])-(4*psi_px[i][j-1]*Z[i][j-1])+(psi_px[i][j-2]*Z[i][j-2])))/(2*DX)-(BETA*psi_px[i][j])+(KAPPA*EDDY[i][j]);
            else if(i==Nx-1 && j<Ny-1)     //forward
                zeta_pt[i][j]=-(((-3*psi_py[i][j]*Z[i][j])+(4*psi_py[i+1][j]*Z[i+1][j])-(psi_py[i+2][j]*Z[i+2][j]))+((-3*psi_px[i][j]*Z[i][j])+(4*psi_px[i][j+1]*Z[i][j+1])-(psi_px[i][j+2]*Z[i][j+2])))/(2*DX)-(BETA*psi_px[i][j])+(KAPPA*EDDY[i][j]);
            else if(i==Nx-1 && j==Ny-1)    //backward
                zeta_pt[i][j]=-(((3*psi_py[i][j]*Z[i][j])-(4*psi_py[i-1][j]*Z[i-1][j])+(psi_py[i-2][j]*Z[i-2][j]))+((3*psi_px[i][j]*Z[i][j])-(4*psi_px[i][j-1]*Z[i][j-1])+(psi_px[i][j-2]*Z[i][j-2])))/(2*DX)-(BETA*psi_px[i][j])+(KAPPA*EDDY[i][j]);
            else if(j==Ny-1 && i<Nx-1)     //forward
                zeta_pt[i][j]=-(((-3*psi_py[i][j]*Z[i][j])+(4*psi_py[i+1][j]*Z[i+1][j])-(psi_py[i+2][j]*Z[i+2][j]))+((-3*psi_px[i][j]*Z[i][j])+(4*psi_px[i][j+1]*Z[i][j+1])-(psi_px[i][j+2]*Z[i][j+2])))/(2*DX)-(BETA*psi_px[i][j])+(KAPPA*EDDY[i][j]);
            else                           //centered
                zeta_pt[i][j]=-(((psi_py[i+1][j]*Z[i+1][j])-(psi_py[i-1][j]*Z[i-1][j]))+((psi_px[i][j+1]*Z[i][j+1])-(psi_px[i][j-1]*Z[i][j-1])))/(2*DX)-(BETA*psi_px[i][j])+(KAPPA*EDDY[i][j]);

        }
    }
}

/*
 * To calculate the relative vorticity one step ahead in time, zeta_t is approximated as follows:
 *
 *      (zeta(k+1)-zeta(k))/DT = F(k)
 *   -> zeta(k+1) = zeta(k) + F(k)*DT
 *
 * The k+1 time step of zeta is equal to the current relative velocity plus the value of zeta_t times change in time
 *
 * In int main() the forecast is found by finding everything at time zero then passing the (k+1) rel. vel. into the
 * functions that use zeta in their calculations and repeated until the desired time period is found
 *
 * ie. with time step DT=30 min for a 24-hr forecast, 48 cycles are needed after time zero
 */
void Z_K(double DT)
{
    for(int i=0; i<Nx; i++)
    {
        for(int j=0; j<Ny; j++)
        {
            Z[i][j]=Z[i][j]+(DT*zeta_pt[i][j]);
        }
    }
}

/*
 * print function is used to either print to the terminal or replace cout with ofstream and write out to a file for
 * MATLAB use.
 */
void print_Z ()
{
    //ofstream os;
    //os.open(s);

    for(int i=0; i<Nx; i++)
    {
        for(int j=0; j<Ny; j++)
        {
            cout<<Z[i][j]<<'\t';
        }
        cout<<'\n';
    }
    //os.close();
}
