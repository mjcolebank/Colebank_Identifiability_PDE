/* The sor06.C main program */

// $Id: sor06.C,v 1.13 2010-10-20 15:38:28 mette Exp $

#include "sor06.h"
#include "tools.h"
#include "arteries.h"
#include "junction.h"

extern "C"  void impedance_init_driver_(int *tmstps);

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

// The vessel-network is specified and initialized, and the flow and
// pressures are to be determed. All global constants must be defined
// in the header file. That is the number of dots per cm, and the number
// of vessels.
int main(int argc, char *argv[])
{
    double tstart, tend, finaltime;
    // Declare parameters: MJC 9/2019
    double f1, f2, f3, fs1, fs2, fs3;
    double alpha_L,alpha_R,beta_L,beta_R,lrr_L,lrr_R,rm_L,rm_R;
    double k1r=0, k2r=0, k3r=0;
    int total_vessels, total_terminal, total_conn, num_pts;
    int ST_calc;
    //double ves_angle; // NOTE: Setting this to 1 will use the
    // ST parameters that have been previously generated

    if (argc < 4) {
        fprintf(stdout,"Not enough entries: Only %d inputs. Exiting now.\n",argc);
        return 1;
    }
    ////// Load in fluids parameters
    f1      = atof(argv[1]);
    f2      = atof(argv[2]);
    f3      = atof(argv[3]);
    fs1     = f1;
    fs2     = f2;
    fs3     = f3;
    alpha_L = atof(argv[4]);
    beta_L = atof(argv[5]);
    lrr_L     = atof(argv[6]);
    rm_L      = atof(argv[7]);
    alpha_R = atof(argv[8]);
    beta_R = atof(argv[9]);
    lrr_R     = atof(argv[10]);
    rm_R      = atof(argv[11]);
    int cycles = 1; // NOTE: Can make this a parameter to pass in
    
    
    // Load in network parameters
    total_vessels  = 3;
    total_terminal = 2;
    num_pts        = 10;
    total_conn     = total_vessels-total_terminal;
    nbrves         = total_vessels;
    // A parameter that can alleviate calculating the impedance for each iteration
    ST_calc = 1;//atoi(argv[14]);

   // Define the files that need to be written to
    int fileID = atoi(argv[12]);
    char namepuALL[20];
    sprintf(namepuALL, "output_%d.2d", fileID);
	FILE *fpALL = fopen (namepuALL, "w");

    // Workspace used by bound_bif
    for(int i=0; i<24; i++) fjac[i] = new double[24];

    tstart    = 0.0;            // Starting time.

    // The number of vessels in the network is given when the governing array of
    // vessels is declared.

    impedance_init_driver_(&tmstps);
    Tube   *Arteries[nbrves];                    // Array of blood vessels.


    int conn_rows = 2*total_vessels; // Initialize with double the number of vessels
    int conn_cols = max_D; // This is defined so that we can have trifurcations if needed.
    int connectivity_matrix[conn_rows][conn_cols];
    int terminal_vessels[total_terminal];
        double dimensions_matrix[total_vessels][3];   // Length and radius

    // Begin sequence of loading text files
    FILE *conn;
    conn = fopen("connectivity.txt","rt");
    int parent, daughter1, daughter2, daughter3, r_in;

    // MJC: Define boolean vector for converging flow
    int conv_flag[total_vessels];
    for (int i=0; i<total_vessels; i++)
        conv_flag[i]=0;

    // Check to see if we have the connecitivity file
    int conn_id = 0;
    if (conn == NULL)
    {
        fprintf(stdout,"Error: Connectivity File Does Not Exist \n");
        return 1;
    }
    while ((r_in = fscanf(conn, "%d %d %d %d", &parent, &daughter1, &daughter2, &daughter3)) != EOF)
    {
        connectivity_matrix[conn_id][0] = parent;
        connectivity_matrix[conn_id][1] = daughter1;
        connectivity_matrix[conn_id][2] = daughter2;
        connectivity_matrix[conn_id][3] = daughter3;
        conn_id++;

        if (daughter1<0)
        {
            fprintf(stdout, "negative daughter 1 %d\n", daughter1);
            conv_flag[-1*daughter1] = -1;
        }
        if (daughter2<0)
        {
            fprintf(stdout, "negative daughter 2 %d\n", daughter2);
            conv_flag[-1*daughter2] = -1;
        }
        if (daughter3<0)
        {
            fprintf(stdout, "negative daughter 3 %d\n", daughter3);
            conv_flag[-1*daughter3] = -1;
        }
    }
    fclose(conn);
    size_conn = conn_id-1;


    dimensions_matrix[0][0] = 4.30;
    dimensions_matrix[1][0] = 2.50;
    dimensions_matrix[2][0] = 5.75;

    dimensions_matrix[0][1] = 1.35;
    dimensions_matrix[1][1] = 0.90;
    dimensions_matrix[2][1] = 1.10;

    dimensions_matrix[0][2] = 1.35;
    dimensions_matrix[1][2] = 0.90;
    dimensions_matrix[2][2] = 1.10;

    terminal_vessels[0] = 1;
    terminal_vessels[1] = 2;



    /* Initialization of the Arteries.
    * // The tube class takes in the following values (in this order)
    * Length: the length of the artery
    *
    * Top radius: the radius at the entry of the vessel
    *
    * Bottom radius: the radius at the end of the vessel
    *
    * Daughter Pointer: the pointer to the connectivity of all the blood
    * vessels, passed to junction.c
    *
    * Minimum radius: minimum radius (in ST)
    *
    * Num Points: the number of pts (per non-dimensional length) to use in numerical solution
    *
    * Init: Set to 1 to make this vessel the inlet vessel (i.e., it gets a flow prescribed); set 0 else
    *
    * f1,f2,f3: Stiffness in the large blood vessels (Eh/r0 = f1*exp(f2*r0)+f3)
    *
    * fs1,fs2,fs3: Stiffness in the small blood vessels (Eh/r0 = fs1*exp(fs2*r0)+fs3) (in ST)
    *
    * STpar_1: the first parameter describing microvascular branching (could be alpha or asym) (in ST)
    *
    * STpar_2: the second parameter describing microvascular branching (could be beta or expo) (in ST)
    *
    * Length-radius ratio: the third parameter for microvasculature; relates length to radius (in ST)
    *
    * Terminal-ID: ID of the terminal vessel, which is utilized when simultations are run without
    * recalculating the structured tree parameters (in ST)
    *
    * ST_calc: A flag to tell the code whether to recalculate the ST parameters (set to 0 to not calculate ST value)
    */

    // MJC: For trying to pass junction condition around
    int *daughter_ptr = &connectivity_matrix[0][0];
    int term_id = total_terminal-1;
    conn_id = total_conn-1;//+total_sten;

    

    // Fixed network geometry
    Arteries[2] = new Tube( dimensions_matrix[2][0], dimensions_matrix[2][1], dimensions_matrix[2][2], k1r, k2r, k3r,
                                               daughter_ptr, rm_L, num_pts, 0, f1,f2,f3,fs1,fs2,fs3, alpha_L, beta_L, lrr_L,term_id,ST_calc, 90);
    
    Arteries[1] = new Tube( dimensions_matrix[1][0], dimensions_matrix[1][1], dimensions_matrix[1][2], k1r, k2r, k3r,
                                               daughter_ptr, rm_R, num_pts, 0, f1,f2,f3,fs1,fs2,fs3, alpha_R, beta_R, lrr_R,term_id,ST_calc, 90);
   
    Arteries[0] = new Tube( dimensions_matrix[0][0], dimensions_matrix[0][1], dimensions_matrix[0][2], k1r, k2r, k3r,
                                       daughter_ptr, 0, num_pts, 1, f1,f2,f3,fs1,fs2,fs3, 0, 0, 0,0,ST_calc, 90);
   
    
    // Solves the equations until time equals tend.///////////////////////////
    /* ADDED BY MJ Colebank
    * Rather than specifying the number of cycles as an input to the function,
    * we want to test to see if the solution has converged. If so, we should exit.*/

    int period_counter = 1; // Count the number of periods you have solved for
    double norm_sol = 1e+6;
    double sol_tol  = 1e-2;
//     fprintf(stdout,"NORM_SOL: %f\n",norm_sol);
    double sol_p1[tmstps],sol_p2[tmstps];
    tend      = Deltat;


    // SOLVE THE MODEL ONCE
    // Note: Only want to test the pressure at the inlet
    int sol_ID = 0;
    while (tend<=period_counter*Period)
    {
        solver (Arteries, tstart, tend, k, WALL_MODEL);
        sol_p1[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0],WALL_MODEL); // for printing
        sol_p1[sol_ID] *= rho*g*Lr/conv;
        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
        sol_ID++;
    }


    // LOOP FOR CONVERGENCE
    double sse;
    while (norm_sol>=sol_tol)
    {
        sol_ID = 0;
        sse    = 0;
        period_counter++;
        while (tend<=period_counter*Period && period_counter<=50)
        {
            solver (Arteries, tstart, tend, k, WALL_MODEL);
            sol_p2[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0],WALL_MODEL); // for printing
            sol_p2[sol_ID] *= rho*g*Lr/conv;
            sse = sse+ sq(sol_p1[sol_ID]-sol_p2[sol_ID]);
            tstart = tend;
            tend   = tend + Deltat; // The current ending time is increased by Deltat.
            sol_ID++;
        }
        norm_sol = sse;
        memcpy (sol_p1, sol_p2, sizeof(sol_p2));
//         printf("NORM_SOL:%f\n",norm_sol);
        fflush(stdout);
    }
//     printf("num_cylces:%d\n",period_counter);
    fflush(stdout);


    // The loop is continued until the final time
    // is reached. If one wants to make a plot of
    // the solution versus x, tend is set to final-
    // time in the above declaration.

    period_counter++;
    finaltime = (period_counter+(cycles-1))*Period;
    while (tend <= finaltime)
    {
        for (int j=0; j<nbrves; j++)
        {
            int ArtjN = Arteries[j]->N;
            for (int i=0; i<ArtjN; i++)
            {
                Arteries[j]->Qprv[i+1] = Arteries[j]->Qnew[i+1];
                Arteries[j]->Aprv[i+1] = Arteries[j]->Anew[i+1];
            }
        }

        // Solves the equations until time equals tend.
        solver (Arteries, tstart, tend, k, WALL_MODEL);


        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (fpALL, tend, 0,WALL_MODEL);
        }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (fpALL, tend, Arteries[save_id]->N/2,WALL_MODEL);
        }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (fpALL, tend, Arteries[save_id]->N,WALL_MODEL);
        }
        // The time within each print is increased.
        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
        fflush(fpALL);
    }



    // In order to termate the program correctly the vessel network and hence
    // all the vessels and their workspace are deleted.
    for (int i=0; i<nbrves; i++) delete Arteries[i];

    // Matrices and arrays are deleted
    for (int i=0; i<18; i++) delete[] fjac[i];

    //    fclose(fp1);
    //    fclose(fp2);
    //    fclose(fp3);
    //    fclose(fp4);
    //    fclose(fp5);
    //    fclose(fp6);
    //    fclose(fp7);
    //    fclose(fp8);
    //    fclose(fp9);
    //    fclose(fp10);
    //    fclose(fp11);
    //    fclose(fp12);
    //    fclose(fp13);
    //    fclose(fp14);
    //    fclose(fp15);
    //    fclose(fp16);
    //    fclose(fp17);
    //    fclose(fp18);
    //    fclose(fp19);
    //    fclose(fp20);
    //    fclose(fp21);
    //    fclose(fp22);
    //    fclose(fp23);
    //    fclose(fp24);
    //    fclose(fp25);
    //    fclose(fp26);
    //    fclose(fp27);
    //    fclose(fp28);
    //    fclose(fp29);
    //    fclose(fp30);
    //    fclose(fp31);
    //    fclose(fp32);
    //    fclose(fp33);
    //    fclose(fp34);
    //    fclose(fp35);
    //    fclose(fp36);
    //    fclose(fp37);
    //    fclose(fp38);
    //    fclose(fp39);
    //    fclose(fp40);
    //    fclose(fp41);
    //    fclose(fp42);
    //    fclose(fp43);
    //    fclose(fp44);
    //    fclose(fp45);
    //    fclose(fp46);
    //    fclose(fp47);
    //    fclose(fp48);
    //    fclose(fp49);
    //    fclose(fp50);
    //    fclose(fp51);
    //    fclose(fp52);
    //    fclose(fp53);
    //    fclose(fp54);
    //    fclose(fp55);
    //    fclose(fp56);

    fclose(fpALL);

    return 0;
}
