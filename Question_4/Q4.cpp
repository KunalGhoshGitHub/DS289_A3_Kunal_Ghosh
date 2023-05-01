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

// This function calculates the analytical solution
void analytical_solution(vector<real> &x, real P, real L, real Ac, real E, vector <real> &Analytical_solution)
{
    for (int i = 0;i < x.size();i++)
    {
        Analytical_solution[i] = (P*x[i]*(x[i]-L))/(2.0*Ac*E);
    }
}

// This function calculate the matrix A of the equation AU = b
void A_matrix(vector<vector<real>> &A, int N)
{
    real rows = N-2;

    A[0][0] = 2.0;
    A[0][1] = -1.0;

    for (int i = 1;i < rows-1;i++)
    {
        A[i][i] = 2.0;
        A[i][i-1] = -1.0;
        A[i][i+1] = -1.0;
    }

    A[rows-1][rows-1] = 2.0;
    A[rows-1][rows-1-1] = -1.0;
}


int main()
{
    vector<real> Value_of_N_s,Value_of_rd_s,Inputs;
    vector<string> Output_file_names;
    real Ac,E,L,P,dx,X_l,X_u,Tolerance,temp;

    // Reading the inputs from the input file
    Inputs = input_parameters("Input.txt");
    Ac = Inputs[0];
    E = Inputs[1];
    L  = Inputs[2];
    P = Inputs[3];
    dx  = Inputs[4];
    Tolerance = Inputs[5];

    // N is the number of the grid points
    int N;
    N = (L/dx)+1.0;

    vector<real> x;
    X_l = 0.0;
    X_u = L;

    // Discretize the domain
    Discretize_1D(x,dx,X_l,X_u);

    // As we have removed the first and last, rows as well as columns
    int rows = N-2;
    int cols = N-2;
    vector<real> u_n(N-2,0.0);
    vector<real> Analytical_solution(N-2,0.0);
    vector<vector <real>> A(rows,vector<real> (cols,0.0));

    // Calculating the A matrix
    A_matrix(A,N);

    // If A matrix is required to be written on a file then uncomment the line below
    // write_to_file(A,"A_matrix.csv");

    // Calculating the RHS
    temp = -P*dx*dx/(Ac*E);

    vector<real> RHS(N-2,temp);

    // If RHS is required to be written on a file then uncomment the line below
    // write_to_file(RHS,"RHS.csv");

    // Using Jacobi Iterative Method to solve the system of equations
    jacobi_solver(A, RHS,Tolerance,u_n);

    // As u_0 = 0 (Boundary Condition)
    u_n.insert(u_n.begin(),0.0);
    Analytical_solution.insert(Analytical_solution.begin(),0.0);

    // As u_N = 0 (Boundary Condition)
    u_n.push_back(0.0);
    Analytical_solution.push_back(0.0);

    // Writing the numerical solution to a file
    write_to_file(u_n,"Question_4_u_Numerical_Solution.csv");
    Output_file_names.push_back("Question_4_u_Numerical_Solution.csv");

    // Calculating the analytical solution
    analytical_solution(x,P,L,Ac,E,Analytical_solution);

    // Writing the analytical solution to a file
    write_to_file(u_n,"Question_4_u_Analytical_Solution.csv");
    Output_file_names.push_back("Question_4_u_Analytical_Solution.csv");

    // Writing the names of all the generated files to a file
    write_to_file(Output_file_names,"Output_file_names.csv");
}
