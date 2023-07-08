#include <iostream> 
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "gsl_rng.h"
using namespace std;
const int N = 100; // Number of particles
const int M = 1; // Number of reference images used

gsl_rng *tau;

// Adrián Marín Boyero, Física UGR

//---------------------------------------------------------------------------------------------------
/*

This C++ script processes storages patterns and storages it into the Hopfield algorithm.

For more than one reference images, they are all storaged in a 2D-array named "references"

To change the functionality of the script, you can change the string "modo" down below. 

*/
//---------------------------------------------------------------------------------------------------

class image
{public:
    int data[N][N];
    ifstream input;
    string name;
};

void generate_s(int s[N][N], gsl_rng *tau, string mode) // We generate our initial set of 1s and 0s
{
    double sesgo = -0.5;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            if(mode == "sesgo")
            { 
                if( (gsl_rng_uniform(tau) - sesgo) > 0.5)
                {
                  s[i][j] = 0;
                }
                else{ s[i][j] = 1;}                
            }

            if(mode == "random")
            {
                s[i][j] = gsl_rng_uniform_int (tau,2);
            }

            if(mode == "1s"){ s[i][j] = 1; }
            if(mode == "0s"){ s[i][j] = 0; }
        }}
}

void print_data_for_python(int data[][N], ofstream& output, string file_name, int MS) // We print an N*N array to be read by the Python script
{
    output << MS << endl;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            output << data[i][j];
            if(j<N-1)
            {
                output << ", ";
            }} 
        output << endl;        
    }
    output << endl;
}

void get_data_for_a_binary_image(int data[][N], ifstream& input, string file_name) // We extract data from reference binary images
{
    input.open(file_name);
    if (input.is_open()){cout << "Succesfully opened: " + file_name << endl;}
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            input >> data[i][j];
        }}
    input.close();
}

 void combine_data_into_one_array(int references[N][N*M], int data_1[][N], int data_2[][N], int data_3[][N], int data_4[][N]) // We combine al ref. images' data into one 2d-rray
 {
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            references[i][j] = data_1[i][j]; //*
            if(M > 1)
            {
            references[i][j+N] = data_2[i][j]; 
            }
            if(M > 2)
            {
            references[i][j+N*2] = data_3[i][j]; 
            }
            if(M > 3)
            {
            references[i][j+N*3] = data_4[i][j]; 
            }
            
        }}}

double calculate_alpha(int references[N][N*M], int m) // We calcualte alpha a specific ref. image's data
{
    double alpha = 0;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            alpha += (double) references[i][j+N*(m)];
        }}
    return 1.0/(double(N*N)) * alpha;
}

double calculate_omega(int i, int j, int k, int l, int references[N][N*M], double alpha[M]) // We calculate synaptic weights (ω)
{
    double ω = 0;
    if(!((i == k) && (j == l)))
    {
        for(int m=0; m<M; m++)
        {
            ω += ((references[i][j + m*N] - alpha[m]) * (references[k][l + m*N] - alpha[m]));
        }}
    return ω * 1/(double(N*N));
}

double calculate_theta(int references[N][N*M], int i, int j, double alpha[M+1]) // We calculate the trigger threshold (θ)
{
    double theta = 0; 
    for(int k=0; k<N; k++) 
    {
        for(int l=0; l<N; l++)
        {
            theta += calculate_omega(i, j, k, l, references, alpha);
        }
    }
    return 0.5 * theta;
}

double calculate_energy_difference(int s[][N], int references[N][N*M], int i, int j, double alpha[M+1]) // Calculation of the diff. in energy
{
    double energy_difference, Δs;
    if(s[i][j] == 1) // We define the difference in value for Δs
    { Δs = -1; }
    else
    { Δs = +1; }
    double sum = 0;
    for(int n=0; n<N; n++) // We calculate the sums of ω_{ij, kl} * s[k][l]
    {
        for(int m=0; m<N; m++) 
        {
            sum += calculate_omega( i, j, n, m, references, alpha) * (double) s[n][m];
        }
    }
    energy_difference = - 2 * Δs * sum  - 2 * calculate_theta(references, i, j, alpha)*Δs;
    return energy_difference;
}

void montecarlo(int s[N][N], int references[N][N*M], double alpha[M+1], gsl_rng *tau, double T) // Montecarlo method
{
    int i = int(1000 * gsl_rng_uniform(tau)) % N, j = int(1000 * gsl_rng_uniform(tau))%N;
    double p, energy_difference;

    energy_difference = calculate_energy_difference(s, references, i, j, alpha);

    if(exp(-energy_difference/T) > 1)
    { p = 1;}
    else
    { p = exp(-energy_difference/T); }

    double random_number = gsl_rng_uniform(tau);
    if(random_number < p)
    {
        s[i][j] = 1 - s[i][j];
    }
}

void overlapping(int s[][N], int references[N][N*M], double alpha[M], int MS, double& overlap, int m)
{
    for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                overlap += (1 / (pow(N, 2) * alpha[m] * (1- alpha[m]) ) ) * (references[i][j + m*N] - alpha[m]) * (s[i][j] - alpha[m]);
            }
        }
}

void calculate_and_print_overlapping(int s[][N], int references[N][N*M], double alpha[M], int MS, ofstream& out_overlapping) // Calculation and printing of overlapping
{
    double overlap = 0;
    cout << " " << fixed << setprecision(4) << MS << " : ";
    out_overlapping << MS << " ";
    for(int m=0; m<M; m++)
    {
        overlapping(s, references, alpha, MS, overlap, m);
        cout << overlap << "  ";
        out_overlapping << fixed << setprecision(6) << overlap << " "; 
        overlap = 0;
    }
    out_overlapping << endl;
    cout << endl;
}

void repetitions_montecarlo(int s[][N], int references[N][N*M], ofstream& output, string output_name, gsl_rng *tau, int rep, string modo, double T, ofstream& out_T) // Sequence of montecarlo's steps
{
    double alpha[M];
    for(int i=0; i<M; i++)
    {
        alpha[i] = calculate_alpha(references, i);
    }
    
    ofstream out_overlapping;

    int total_steps = rep*N*N, show = 0, MS = 0;

    output.open(output_name);
    out_overlapping.open("overlapping.txt");
    
    if(modo == "normal")
    {
        cout << endl << "Overlapping (m) -----------------------" << endl << " MS   image 1  &  image 2  &  image 3  &  image 4:" << endl;
        for(int steps=0; steps<=total_steps; steps++)
        {
            montecarlo(s, references, alpha, tau, T);
            if(show == N*N)
            {
                MS++;
                calculate_and_print_overlapping(s, references, alpha, MS, out_overlapping);
                print_data_for_python(s, output, output_name, MS);
                show = 0;
            }
            show++;   
        }

        ofstream final;
        final.open("patron_final.txt");
        print_data_for_python(s, final, "patron_final.txt", MS);
        final.close();
    }

    if(modo == "temperaturas")
    {
        for(int steps=0; steps<=total_steps; steps++)
        {
            montecarlo(s, references, alpha, tau, T); 
        }

        double overlap = 0;
        out_T << fixed << setprecision(5) << T << " ";
        for(int i=0; i<M; i++)
        {
            overlapping(s, references, alpha, MS, overlap, i);
            out_T << fixed << setprecision(5) << overlap << " ";
            overlap = 0;
        }
        out_T << endl;

        ofstream final;
        final.open("patron_final.txt");
        print_data_for_python(s, final, "patron_final.txt", MS);
        final.close();

    }
    
    if(modo == "remembered_patterns")
    {
        for(int steps=0; steps<=total_steps; steps++)
        {
            montecarlo(s, references, alpha, tau, T);
            if(show == N*N)
            {
                MS++;
                print_data_for_python(s, output, output_name, MS);
                show = 0;
            }
            show++;   
        }

        double number_P = 0, overlap = 0;
        for(int i=0; i<M; i++)
        {
            overlapping(s, references, alpha, MS, overlap, i);
            if(abs(overlap) >= 0.75)
            {
                number_P += 1;
            }
            overlap = 0;
        }

        cout << "Para " << M << " patrones se recuerdan " << number_P << endl;

    }

    out_overlapping.close();
    output.close();
}

void deformacion(int s[N][N], int patron[N][N], gsl_rng *tau, double porciento)
{
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            if(gsl_rng_uniform(tau) <= 0.01*porciento)
            {
                s[i][j] = 1 - patron[i][j];
            }
            else
            {
                s[i][j] = patron[i][j];
            }
        }
    }
}

int main()
{
    extern gsl_rng *tau;
    int semilla=18237247;
    tau=gsl_rng_alloc(gsl_rng_taus);

    int s[N][N], references[N][N*M];
    ofstream output; string output_name = "output_hopefield_gsl.txt";
    ofstream imagenes;

    // Number of Montecarlo Steps (MS) we want to run
    int MS = 30;

    string modo = "normal"; //"normal", "temperaturas" o "remembered_patterns" (apartado4)
    double T = pow(10, -4); // Temperature

    if(modo == "normal")
    {

        // We define the reference images we will get binary data from
        image salchicha; salchicha.name = "salchicha.txt";
        image casa_azul; casa_azul.name = "casa_azul.txt";
        image palmera; palmera.name = "palmera.txt";
        image octopus; octopus.name = "octopus.txt";

        // For ortognal images
        /*
        image salchicha; salchicha.name = "salchicha_o.txt";
        image casa_azul; casa_azul.name = "casa_azul_o.txt";
        image palmera; palmera.name = "palmera_o.txt";
        image octopus; octopus.name = "octopus_o.txt";
        */

        // We get data from the ref. images and get them all into a 2d array [N][N*number of images]
        get_data_for_a_binary_image(salchicha.data, salchicha.input, salchicha.name);
        get_data_for_a_binary_image(casa_azul.data, casa_azul.input, casa_azul.name);
        get_data_for_a_binary_image(palmera.data, palmera.input, palmera.name);
        get_data_for_a_binary_image(octopus.data, octopus.input, octopus.name);

        generate_s(s, tau, "sesgo"); // We generate our initial, 2d binary array (random, sesgo vs 1s, 0s)

        //deformacion(s, salchicha.data, tau, 50);

        //imagenes.open("im_salchicha.txt");
        //print_data_for_python(s, imagenes, "im_salchicha.txt", 0);
        //imagenes.close();

        combine_data_into_one_array(references, salchicha.data, casa_azul.data, palmera.data, octopus.data);

        // We do the sequence of montecarlo's steps and print the results
        repetitions_montecarlo(s, references, output, output_name, tau, MS, modo, T, output);

        /*
        imagenes.open("im_salchicha.txt");
        print_data_for_python(salchicha.data, imagenes, "im_salchicha.txt", 0);
        imagenes.close();
        */
 
    }

    if(modo == "temperaturas")
    {
        // We define the reference images we will get binary data from
        image salchicha; salchicha.name = "salchicha.txt";
        image casa_azul; casa_azul.name = "casa_azul.txt";
        image palmera; palmera.name = "palmera.txt";
        image octopus; octopus.name = "octopus.txt";

        // For ortognal images
        /*
        image salchicha; salchicha.name = "salchicha_o.txt";
        image casa_azul; casa_azul.name = "casa_azul_o.txt";
        image palmera; palmera.name = "palmera_o.txt";
        image octopus; octopus.name = "octopus_o.txt";
        */


        // We get data from the ref. images and get them all into a 2d array [N][N*number of images]
        get_data_for_a_binary_image(salchicha.data, salchicha.input, salchicha.name);
        get_data_for_a_binary_image(casa_azul.data, casa_azul.input, casa_azul.name);
        get_data_for_a_binary_image(palmera.data, palmera.input, palmera.name);
        get_data_for_a_binary_image(octopus.data, octopus.input, octopus.name);

        combine_data_into_one_array(references, salchicha.data, casa_azul.data, palmera.data, octopus.data);

        //generate_s(s, tau, "random"); // We generate our initial, 2d binary array (random, sesgo vs 1s, 0s)
        
        deformacion(s, octopus.data, tau, 40);

        ofstream out_T;
        out_T.open("registro_temperaturas.txt");
        imagenes.open("simulacion_temperaturas.txt");

        for(int k=0; k<=3; k++)
        {
            for(int n=2.5; n<=10; n+=2)
            {
                repetitions_montecarlo(s, references, output, output_name, tau, MS, modo, n*T*pow(10, k), out_T);
                print_data_for_python(s, imagenes, "simulacion_temperaturas", n*T*pow(10, k));
                deformacion(s, octopus.data, tau, 40);
            }

        }
        imagenes.close();
        out_T.close();
    
    }

    if(modo == "remembered_patterns")
    {
        int ref_patterns[N][N*M];

        // We extract our random-generated starting condition
        image inicial; inicial.name = "inicial.txt";
        get_data_for_a_binary_image(inicial.data, inicial.input, inicial.name);

        for(int i=0; i<N; i++)
        {
            for(int j=0; j<(N*M); j++)
            {
                ref_patterns[i][j] = gsl_rng_uniform_int(tau, 2);
            }
        }

        repetitions_montecarlo(inicial.data, ref_patterns, output, output_name, tau, MS, modo, T, output);
    }

    return 0;
}