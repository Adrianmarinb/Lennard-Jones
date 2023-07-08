#include <iostream>
#include <fstream>
#include <cmath>
#include "gsl_rng.h"
using namespace std;
gsl_rng *tau;

// --- Constant Parameters -----------------------------------------

const int PartN = 16; // Number of particles
const double sigma = 1; // Parameter sigma
const double ε = 1; // Energy factor
const double h = 0.002; // Time accuracy
const int L = 4; // Length of the simulation box's side

const double rc = 3; // Cutoff distance
const double Fc = (48/pow(rc, 13)) - (24/pow(rc, 7)); // Force felt by two particles at cutoff distance
const double Vc = (4/pow(rc, 12)) - (4/pow(rc, 6)); // Potential felt by two particles at cutoff distance

// --- Introduction -------------------------------------------------

/*

For this code, we will store our data in a 2-dim array, data[PartN][6],
in such a way that we'll have, for a particle i:

data[i][0] = x_position
data[i][1] = y_position
data[i][2] = x_velocity
data[i][3] = y_velocity
data[i][4] = x_acceleration
data[i][5] = y_acceleration

Also, we are assuming that beyond a distance rc between two particles, their interaction is zero. That way, we'll have
a modified force, Fm, which will be the calculated force, F, minus the force of the rc distance, Fc. That way we'll have a continuous function. 

Also, we are using periodic boundary conditions. 

*/

// --- Functions and main sequence ----------------------------------

// Initialization of the data matrix with empty values
void initialize_matrix(long double data[PartN][6])
{
    for (int i = 0; i < PartN; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            data[i][j] = 0;
        }
    }
}

// Function that calculates the closest distance between 2 particles by making use of periodic boundary conditions (LxL box)
long double distance(long double& dx, long double& dy)
{
    if(abs(dx) > L/2)
    {
        dx = (L - abs(dx)) * (-dx)/abs(dx);
    }
    if(abs(dy) > L/2)
    {
        dy = (L - abs(dy)) * (-dy)/abs(dy);
    }

    return sqrt(dx*dx + dy*dy);
}

// Function that generates a random distribution of particles, separated by an arbitrary minimum distance
void random_initial_distribution(long double data[PartN][6], gsl_rng *tau)
{
    double min_distance = 0.7;
    long double dx, dy;
    bool condition = false;

    data[0][0] = gsl_rng_uniform(tau) * L;
    data[0][1] = gsl_rng_uniform(tau) * L;

    for (int i = 1; i < PartN; i++)
    {
        condition = true;
        while (condition == true)
        {
            data[i][0] = gsl_rng_uniform(tau) * L;
            data[i][1] = gsl_rng_uniform(tau) * L;

            condition = false;

            for (int j = i - 1; j >= 0; j--)
            {
                dx = data[i][0] - data[j][0];
                dy = data[i][1] - data[j][1];

                if (distance(dx, dy) < min_distance)
                {
                    condition = true;
                    break;
                }
            }
        }
    }
    cout << "Successful random distribution" << endl; 
}

// Function that generates a aquare-net initial distribution 
void square_initial_distribution(long double data[PartN][6])
{
    int j = 0;
    for (int i = 0; i < PartN; i++)
    {
        if ((i % 4 == 0) && (i != 0))
        {
            j++;
        }
        data[i][0] = L / 8.0 + L / 4.0 * (i % 4);
        data[i][1] = L / 8.0 + L / 4.0 * j;
    }
}

// Function that generates an hexagonal initial distribution 
void honeycomb_initial_distribution(long double data[PartN][6])
{
    int j = 0;
    for (int i = 0; i < PartN; i++)
    {
        if ((i % 4 == 0) && (i != 0))
        {
            j++;
        }
        data[i][0] = L / 8.0 + L / 4.0 * (i % 4) - L / 8.0 * (j % 2);
        data[i][1] = L / 8.0 + L / 4.0 * j;
    }
}

// Function that gives each particle a velocity of a certain module, but in random directions
void random_initial_velocities_module(long double data[PartN][6], gsl_rng *tau, double v_module)
{
    double angle;
    for (int i = 0; i < PartN; i++)
    {
        angle = gsl_rng_uniform(tau) * 2 * 3.141592;
        data[i][2] = cos(angle) * v_module;
        data[i][3] = sin(angle) * v_module; 
       /*
       data[i][2] = gsl_rng_uniform(tau);
       data[i][3] = 0;
       */
    }
}

// Function to print data to be read by the Python script for visualization
void print_data_for_python(long double data[PartN][6], ofstream &output, string file_name, double t)
{
    output << t << endl;
    for (int i = 0; i < PartN; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            output << data[i][j];
            if (j != 1)
            {
                output << ", ";
            }
        }
        output << endl;
    }
    output << endl;
}

// Function that calculates the forces between particles due to the Lennard-Jones potential
void force_calculations(long double data[PartN][6], long double& Potential_Energy)
{
    Potential_Energy = 0;
    for(int k=0; k<PartN; k++)
    {
        for(int l=4; l<6; l++)
        {
            data[k][l] = 0;
        }
    }
    
    long double dx, dy, r, F, Fm, Vm, V;
    for(int i=0; i< (PartN - 1); i++)
    {
        for(int j=i+1; j<PartN; j++)
        {
            dx = data[i][0] - data[j][0];
            dy = data[i][1] - data[j][1];
            r = distance(dx, dy);

            if( r <= rc)
            {
                F = (48/pow(r, 13)) - (24/pow(r, 7));
                Fm = F - Fc;

                data[i][4] += Fm * (dx/r);
                data[i][5] += Fm * (dy/r);
                data[j][4] -= Fm * (dx/r);
                data[j][5] -= Fm * (dy/r);
                V = (4/pow(r, 12)) - (4/pow(r, 6));
                Vm = V - Vc + r*Fc - rc*Fc;
                Potential_Energy += Vm;
            }
        }
    }
}

// Function resposible for updating the position of particles by making use of periodic boundary conditions
void update_positions(long double data[PartN][6], long double omega[PartN][2], double& momentum)
{
    momentum = 0;
    for(int i=0; i<PartN; i++)
    {
        for(int j=0; j<2; j++)
        {
            data[i][j] += h*omega[i][j];

            if(data[i][j] > L)
            {
                momentum += abs(data[i][j + 2]);
                data[i][j] = fmod(data[i][j], L);
            }
            if(data[i][j] < 0)
            {
                momentum += abs(data[i][j + 2]);
                data[i][j] = L - fmod(abs(data[i][j]), L);
            }
        }
    }
}

// Function that updates the velocity of each particle
void update_velocities(long double data[PartN][6], long double omega[PartN][2], long double& Kinetic_Energy)
{
    Kinetic_Energy = 0;
    for(int i=0; i<PartN; i++)
    {
        for(int j=0; j<2; j++)
        {
            data[i][j + 2] = omega[i][j] + (h/2.0) * data[i][j + 4];
        }
        Kinetic_Energy += 0.5*(data[i][2]*data[i][2] + data[i][3]*data[i][3]);
    }
}

// Function that makes use of the velvet method to update de positions, velocities and accelerations of each particle in the system
void velvet_integration(long double data[PartN][6], long double& Potential_Energy,  long double& Kinetic_Energy, double& momentum)
{
    // We define omega[PartN][2] as an auxiliary 2-dim array
    long double omega[PartN][2];

    // We evaluate omega using omega = vel + (h/2)*accel
    for (int i = 0; i < PartN; i++) 
    {
        for (int j = 0; j < 2; j++)
        {
            omega[i][j] = data[i][j + 2] + (h / 2) * data[i][j + 4];
        }
    }

    update_positions(data, omega, momentum);
    force_calculations(data, Potential_Energy);
    update_velocities(data, omega, Kinetic_Energy);
}

// This function serves the purpose of treating data related to some parameters of the simulation, as well as printing them onto a .txt
void print_other_parameters(long double data[PartN][6], double t, long double Potential_E, long double Kinetic_E, ofstream &output, double& momentum, double initial_r_0[2])
{
    long double Average_PE, Average_KE, Total_E, sum_of_v = 0, fluctuation, dx, dy, distance_1p, distance_2p, distance_3p;

    // We calculate the average PE, KE and Total Energy
    Average_PE = Potential_E/PartN;
    Average_KE = Kinetic_E/PartN;
    Total_E = Average_PE + Average_KE;

    // We calculate the distance between particle 1 and its initial position
    dx = data[0][0] - initial_r_0[0];
    dy = data[0][1] - initial_r_0[1];
    fluctuation = pow(distance(dx, dy), 2);

    // We calculate the distance between particles 1 and 2
    dx = data[0][0] - data[1][0];
    dy = data[0][1] - data[1][1];
    distance_1p = pow(distance(dx, dy), 2);

    dx = data[0][0] - data[5][0];
    dy = data[0][1] - data[5][1];
    distance_2p = pow(distance(dx, dy), 2);

    dx = data[0][0] - data[7][0];
    dy = data[0][1] - data[7][1];
    distance_3p = pow(distance(dx, dy), 2);

    // We calculate the mean velocity of the particles
    for(int i=0; i<PartN; i++)
    {
        sum_of_v += pow(data[i][2], 2) + pow(data[i][3], 2);
    }
    sum_of_v *= 1.0/PartN;

    // We print the data for the Python script to read
    output << t << " " << Average_PE << " " << Average_KE << " " << Total_E << " " << sum_of_v << " " << 2.0 * momentum  << " " << fluctuation << " " << distance_1p <<  " " << distance_2p << " " << distance_3p << endl;
}

// This function registers the |V|, Vx, and Vy of the particles in a series of intervals, so as to be represented in bars by the Python script
void v_statistics(long double data[PartN][6], double v_distribution[100][4])
{
    double v_module;
    double size = 0.2; // Size of the intervals

    for(int k=0; k<100; k++)
    {
        v_distribution[k][0] = (k-50)*size;
    }

    // We register the velocities of each particle (|V|, Vx and Vy)
    for(int i=0; i<PartN; i++)
    {
        v_module = sqrt(pow(data[i][2], 2) + pow(data[i][3], 2));
        
        for(int j=0; j<100; j++)
        {
            if( (v_module >= size*(j-50)) && (v_module < size*(j-49)))
            {
                v_distribution[j][1] += 1;
            }
            if( (data[i][2] >= size*(j-50)) && (data[i][2] < size*(j-49)))
            {
                v_distribution[j][2] += 1;
            }
            if( (data[i][3] >= size*(j-50)) && (data[i][3] < size*(j-49)))
            {
                v_distribution[j][3] += 1;
            }
        }
    }
}

// This function serves the purpose of treating data related to particles' velocities, as well as printing it onto a .txt
void print_v_stats(double v_distribution[100][4])
{
    ofstream v_stats; 
    double sum_v = 0, sum_vx = 0, sum_vy = 0;

    for(int i=0; i<100; i++)
    {
        sum_v += v_distribution[i][1];
        sum_vx += v_distribution[i][2];
        sum_vy += v_distribution[i][3];
    }
    
    // We normalize the total data collected by taken into account the total number of values and the size of the intervals (0.2)
    for(int j=0; j<100; j++)
    {
        v_distribution[j][1] *= 1.0/(sum_v * 0.2) ;
        v_distribution[j][2] *= 1.0/(sum_vx * 0.2);
        v_distribution[j][3] *= 1.0/(sum_vy * 0.2);
    }

    cout << "Sum_v = " << sum_v << " | Sum_vx = " << sum_vx << " | Sum_vy = " << sum_vy << endl;

    // We print the data for the Python script to read
    v_stats.open("v_statistics.txt");
    for(int k=0; k<100; k++)
    {
        v_stats << v_distribution[k][0] << " " << v_distribution[k][1] << " " << v_distribution[k][2] << " " << v_distribution[k][3] << endl;
    }
    v_stats.close();

}

// This function increases the velocity of the particles by a certain factor, at specific timestamps 
void increase_velocity(long double data[PartN][6], double t)
{
    double factor = 1.1;

    bool condition = false;
    for(int j=0; j<20; j++)
    {
        if(condition == false)
        {
            if( abs(t - (20 + j*5)) < h/2.0 )
            //if( ((t >= 20) && (t < 20.002)) || ((t >= 30) && (t < 30.002)) || ((t >= 35) && (t < 35.002)) || ((t >= 45) && (t < 45.002)))
            {
                condition = true;
                for(int i=0; i<PartN; i++)
                {
                    data[i][2] *= factor;
                    data[i][3] *= factor;
                }

                cout << "Velocity changed at t: " << t << " by a factor of " << factor << endl;
            }
        }
    }
}

// This function is the pair correlation function, which studies the density of particles per unit distance
void g_function(long double data[PartN][6], double g[60][2])
{
    double initial_r = 1, final_r = 3.5, n_bins = 60, distance_2p;
    double size = (final_r - initial_r) / n_bins;
    long double dx, dy;

    for(int k=0; k<n_bins; k++)
    {
        g[k][0] = initial_r + k*size;
    }

    for(int i=0; i<PartN; i++)
    {
        for(int j=0; j<PartN; j++)
        {
            dx = data[i][0] - data[j][0];
            dy = data[i][1] - data[j][1];
            distance_2p = pow(distance(dx, dy), 2);
            for(int m=0; m<n_bins; m++)
            {
                if( (distance_2p >= (initial_r + m*size)) && (distance_2p < (initial_r + (m+1)*size)) )
                {
                    g[m][1] +=1;  
                }              
            }
        }
    }
}

// This function serves the purpose of treating data related to the g_function as well as printing it onto a .txt
void print_g_stats(double g[60][2])
{
    ofstream g_stats;
    double sum_g = 0;

    for(int i=0; i<60; i++)
    {
        sum_g += g[i][1];
    }

    // We normalize the total data collected by taken into account the total number of values and the size of the intervals (0.1)
    for(int j=0; j<60; j++)
    {
        g[j][1] *= 1/(sum_g * 0.1);
    }

    // We print the data for the Python script to read
    g_stats.open("g_statistics.txt");
    for(int k=0; k<60; k++)
    {
        g_stats << g[k][0] << " " << g[k][1] << endl;
    }
    g_stats.close();

}

// General algorithm for the iterations carried out for the simulations
void general_algorithm(long double data[PartN][6], ofstream &output, string output_name, int total_iterations, double initial_r_0[2])
{
    ofstream parameters_output;
    double t = 0.0, v_distribution[100][4], momentum, g[60][2]; int show = 0;

    // We initilalize the v_distribution array (for future velocity analysis)
    for(int k=0; k<100; k++)
    {
        v_distribution[k][1] = 0;
        v_distribution[k][2] = 0;
        v_distribution[k][3] = 0;
    }
    // We initilalize the g array (for the pair correlation function)
    for(int i=0; i<60; i++)
    {
        g[i][1] = 0;
    }

    long double Potential_Energy, Kinetic_Energy;

    output.open(output_name);
    parameters_output.open("parameters.txt");

    // We print our initial data into the .txt
    print_data_for_python(data, output, output_name, t);
    force_calculations(data, Potential_Energy);

    // We carry out the following functions a number of arbitrary iterations
    for (int steps = 0; steps < total_iterations; steps++)
    {
        if((t > 20) && (t < 50))
        {
            v_statistics(data, v_distribution);
        }
        velvet_integration(data, Potential_Energy, Kinetic_Energy, momentum);
        print_other_parameters(data, t, Potential_Energy, Kinetic_Energy, parameters_output, momentum, initial_r_0);

        // Interval time at which we want to study the velocities of the particles
        //if((t > 20) && (t < 50))

        // Interval time at which we want to study the pair correlation function
        if((t > 1) && (t < 30))
        {
            g_function(data, g);
        }

        // We increase the time value
        t += h;

        // We print the state of the system each x iterations (so as to avoid an excesively long .txt)
        if (show == 20)
        {
            print_data_for_python(data, output, output_name, t);
            show = 0;
        }
        show++;

        //We update the velocities of the particles at certain timestamps to study some of the system's properties
        increase_velocity(data, t);
    }

    parameters_output.close();
    output.close();

    cout << "The video shows an evolution of " << t << " seconds." << endl;

    // We print the data for particles' velocites and the pair correlation functions
    print_v_stats(v_distribution);
    print_g_stats(g);
    
}

int main()
{
    ofstream output;
    string output_name = "lennard-jones-output.txt";

    extern gsl_rng *tau;
    int seed = 18237247;
    tau = gsl_rng_alloc(gsl_rng_taus);

    long double data[PartN][6];
    initialize_matrix(data);

    // --- Possible initial configurations -----------------------

    //random_initial_distribution(data, tau);
    square_initial_distribution(data);
    //honeycomb_initial_distribution(data);

    double v_module = 0; // Module of the velocities of each particle
    random_initial_velocities_module(data, tau, v_module);

    double initial_r_0[2];
    initial_r_0[0] = data[0][0]; initial_r_0[1] = data[0][1];
    
    // --- Iterations and general algorithm -----------------------

    int total_iterations = 50000;
    general_algorithm(data, output, output_name, total_iterations, initial_r_0); 

    return 0;
}