/* Método de solución de ecuaciones diferenciales ordinarias
    -Runge Kuta Generalizado para sistema de m ecuaciones 
    
    Métodos Numéricos
    Edgar Aguilera Hernandez
    11/23/2024
*/

#include "tools.h"
#include "matrixsltn.h"
#include "factrmethods.h"


void plotGNU(int total_size,double *range){
    //iniciar aplicacion GNUPLOT
    FILE *fPlotter = popen("gnuplot -persist", "w");
    
    //Configurar plot
    fprintf(fPlotter,"set terminal x11 0 size 1200,600\n");
    fprintf(fPlotter,"set title \"Potencial de acción\" \n");
    fprintf(fPlotter,"set xlabel \" tiempo(s) \" \n");
    fprintf(fPlotter,"set xrange [%f:%f] \n",range[0],range[1]);
    fprintf(fPlotter,"set ylabel \" V(mv) \" \n");
    fprintf(fPlotter,"set yrange [-4:4] \n");

    //crear grafica
    fprintf(fPlotter, "plot \"y_sol.txt\" using 1:2 smooth csplines lw 3 lt 6 lc rgb \"blue\" title \"Nodo Sinoatrial \" \n");
    //fprintf(fPlotter, "\"l_sol.txt\" using 1:2 smooth csplines lw 3 lt 4 lc rgb \"red\" title \"Sol. Linces.\",");
    //fprintf(fPlotter, "\"c_data.txt\" pt 6 ps 1.5 lc rgb \"blue\" title \"Datos Conejos.\",");
    //fprintf(fPlotter, "\"l_data.txt\" pt 7 ps 1.5 lc rgb \"red\" title \"Datos Linces.\"\n");
    fprintf(fPlotter,"set terminal png size 1200,600\n");
    fprintf(fPlotter,"set output \"vdp_curve%d.png\"\n",total_size);
    fprintf(fPlotter,"replot\n");

    //phase diagram
    fprintf(fPlotter,"set terminal x11 1 size 1200,600\n");
    fprintf(fPlotter,"set title \"Plano de fase\" \n");
    fprintf(fPlotter,"set xlabel \" y(mv) \" \n");
    fprintf(fPlotter,"set xrange [-4:4] \n");
    fprintf(fPlotter,"set ylabel \" v(t) \" \n");
    //fprintf(fPlotter,"set yrange [-10:10] \n");

    fprintf(fPlotter, "plot \"phase_data.txt\" with lines lw 3 lt 2 title \"Orbita de Sol.\"\n");
    fprintf(fPlotter,"set terminal png size 1200,600\n");
    fprintf(fPlotter,"set output \"phase_diagram%d.png\"\n",total_size);
    fprintf(fPlotter,"replot\n");
    
    fclose(fPlotter);
}

/*Van Der Pol equation as 2 first grade equations with variables t,y,u*/
double vdpOscillator(double *var, int equation){
    //parameters
    double k = 30.0;     //damping factor 
    double a = 60.0;     //restoring factor
    double solution;

    if(equation == 0){
        // y'= u
        solution = var[1];
    }else if(equation == 1){
        //u' = k(1-y^2)u - ay
        solution = k * (1 - var[0]*var[0])*var[1] - a*var[0];
    }

    return solution;
}

/*Van Der Pol modified equation as 2 first grade euqations with variables t,y,u*/
double vdpOscillatorModified(double *var, int equation){
    //parameters
    double v[2] = {1.5,-1.5};   //depolarization - fractal speed 
    double k = 4.0;     //damping factor 
    double d = 6.0;
    double e = 8.0; 
    double solution;

    if(equation == 0){
        // y'= v
        solution = var[1];
    }else if(equation == 1){
        //v' = -k(y - v_1)(y  - v_2)u - [y(y + d)(y + e)] / de
        solution = -k * (var[0]-v[0])*(var[0]-v[1] )*var[1] - 
        (var[0]*(var[0] + d)*(var[0] + e)) / d*e;
    }

    return solution;
}

/*4th order Runge-Kuta method for m ordinary differential equations*/
double **solveRungeKuta(double (*F)(double *,int), double *range, int n, int m, double *init_c){
    //solutions matrix (n+1 solutions, m variables + independent variable)
    double **solution = genMatriz_double(n+1,m+1);
    //set of 4th grade approximations, m equations
    double k[4][m];
    //set size of steps
    double h = (range[1] - range[0]) /n;
    //set starting value 
    double t = range[0];
    solution[0][0] = t;
    //set initial solution as the initial condition
    double w[m];
    double w_aux[m];
    copyVector(init_c,w,m);
    for(int i = 0; i< m; i++){
        solution[0][i+1] = w[i];
    }
    

    //iterate solution approximation over the n steps of [a,b]
    for(int i = 0; i < n; i++){
        //start series of averages to reach 4th-order
        //k-1
        for(int j = 0; j < m; j++){
            k[0][j] = h*(*F)(w,j);
        }
        //add previous estimations to the variables (excluding indepedent variable)
        for(int l = 0; l < m; l++){
            w_aux[l] = w[l] + (0.5)*k[0][l]; 
        } 
        //k-2
        for(int j = 0; j < m; j++){
            k[1][j] = h*(*F)(w_aux,j);
        }
        //add previous estimations to the variables (excluding indepedent variable)
        for(int l = 0; l < m; l++){
            w_aux[l] = w[l] + (0.5)*k[1][l]; 
        } 
        //k-3
        for(int j = 0; j < m; j++){
            k[2][j] = h*(*F)(w_aux,j);
        }
        //add previous estimations to the variables (excluding indepedent variable)
        for(int l = 0; l < m; l++){
            w_aux[l] = w[l] + k[2][l]; 
        } 
        //k-4
        for(int j = 0; j < m; j++){
            k[3][j] = h*(*F)(w_aux,j);
        }

        //update solution of w with k-th estimations
        for(int j = 0; j < m; j++){
            w[j] += (k[0][j] + 2*k[1][j] + 2*k[2][j] + k[3][j])/6.0;
        }
        
        //update t with next step
        t += h; 

        //store solution set of variables for t value
        solution[i+1][0] = t;
        for(int j = 0; j< m; j++){
            solution[i+1][j+1] = w[j];
        }
    }

    return solution;
}

int main(){
    //equation system parameters
    double range[2] = {0,3};   //solution range 
    int n = 500;                 //number of steps 
    int m = 2;                  //number of equations/variables
    double init_c[2] = {2.0,0.0};  //initial conditions y(0), v(0)

    //solve ordinal differential equation
    double **solution = solveRungeKuta(vdpOscillatorModified,range,n,m,init_c);

    //write solution into file to plot
    //write data points to plot
    FILE *fData_y = fopen("y_sol.txt","w");
    FILE *fData_p = fopen("phase_data.txt","w");
    for(int i = 0; i < n+1; i++){
        fprintf(fData_y,"%0.4f %0.4f\n",solution[i][0],solution[i][1]);
        fprintf(fData_p,"%0.4f %0.4f\n",solution[i][1],solution[i][2]);
    }
    fclose(fData_y);
    fclose(fData_p);

    //plot interpolated data
    plotGNU(n,range);

    freeMatrix(1,solution);
}