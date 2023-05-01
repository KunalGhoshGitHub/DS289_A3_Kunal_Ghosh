#include <iostream>
#include <cassert>
#include<algorithm>
#include<vector>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include<functional>
#include "DS289.h"

using namespace std;

typedef double real;

int main()
{
    vector<real> Value_of_N_s,Value_of_rd_s,Inputs;
    vector<string> Output_file_names;
    real pi,dx,X_l,X_u,Tolerance,dt,T_end,Alpha,m,p;

    pi = 4.0*atan(1.0);

    // Reading the inputs from the input file
    Inputs = input_parameters("Input.txt");
    X_l = Inputs[0];
    X_u = Inputs[1];
    Alpha = Inputs[2];
    dt = Inputs[3];
    T_end = Inputs[4];

    // Reading the values of N from the file
    Value_of_N_s = input_parameters("Value_of_Ns.txt");

    // This loop iterates over the different values of N
    for (int n = 0;n < Value_of_N_s.size();n++)
    {
        int N = Value_of_N_s[n];

        vector<real> x;

        // Discretize the domain
        dx = (X_u-X_l)/(N-1.0);
        Discretize_1D(x,dx,X_l,X_u);

        m = dt/(2.0*dx);
        p = Alpha*dt/(dx*dx);

        // Calculating the number of the time steps required
        int num_steps = T_end/(dt);

        vector<real> u_n;
        vector<real> u_n_1(N,0.0);

        // Adjustments for the last time step
        real dt_new = T_end - (dt*num_steps);
        real m_new = dt_new/(2.0*dx);
        real p_new = Alpha*dt_new/(dx*dx);

        // Initial Conditions
        for (int i = 0;i < x.size()-1;i++)
        {
            u_n.push_back(sin(4.0*pi*x[i]) + sin(6.0*pi*x[i]) + sin(10.0*pi*x[i]));
        }

        // Time marching loop
        for (int t_steps = 0; t_steps < num_steps; t_steps++ )
        {
            for (int j = 1;j < N-1;j++)
            {

                u_n_1[j] = u_n[j] - (m*u_n[j]*(u_n[j+1] - u_n[j-1])) + (p*(u_n[j+1] - (2.0*u_n[j]) + u_n[j-1]));
            }
            // Periodic Boundary Conditions
            u_n_1[0] = u_n[0];
            u_n_1[N-1] = u_n[N-1];
            u_n = u_n_1;

            // Storing the solution of some of the intermediate time steps
            if (t_steps%30 == 0)
            {
                write_to_file(u_n_1,"Question_1_b_u_n_t_"+to_string(t_steps*dt)+"_N_"+to_string(N)+"_.csv");
                Output_file_names.push_back("Question_1_b_u_n_t_"+to_string(t_steps*dt)+"_N_"+to_string(N)+"_.csv");
            }
        }

        // Last time step
        for (int j = 1;j < N-1;j++)
        {
            u_n_1[j] = u_n[j] - (m_new*u_n[j]*(u_n[j+1] - u_n[j-1])) + (p_new*(u_n[j+1] - (2.0*u_n[j]) + u_n[j-1]));
        }

        // Periodic Boundary Conditions
        u_n_1[0] = u_n[0];
        u_n_1[N-1] = u_n[N-1];

        // Writing the numerical solution to a file
        write_to_file(u_n_1,"Question_1_b_u_n_t_"+to_string(T_end)+"_N_"+to_string(N)+"_.csv");
        Output_file_names.push_back("Question_1_b_u_n_t_"+to_string(T_end)+"_N_"+to_string(N)+"_.csv");
    }

    // Writing the names of all the generated files to a file
    write_to_file(Output_file_names,"Output_file_names.csv");
}
