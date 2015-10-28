//
//  koptim.cpp
//  ktools
//
//  Created by Joowhi Lee on 9/3/15.
//
//

#include "koptim.h"

#include <iostream>
#include <math.h>
#include <nlopt.h>

using namespace std;

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

typedef struct {
    double a, b;
} my_constraint_data;

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

int test1() {
    double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
    nlopt_opt opt;
    
    opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_min_objective(opt, myfunc, NULL);
    
    my_constraint_data data[2] = { {2,0}, {-1,1} };

    nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
    nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);

    nlopt_set_xtol_rel(opt, 1e-4);

    double x[2] = { 1.234, 5.678 };  /* some initial guess */
    double minf; /* the minimum objective value, upon return */
    
    if (nlopt_optimize(opt, x, &minf) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
    }
    
    nlopt_destroy(opt);
	return 0;
}



double hfunc(unsigned n, const double *x, double *grad, void *my_func_data) {
	static int count = 0;

	count ++;
	if (count % 1000 == 0) {
		cout << "iteration: " << count << endl;
	}
    

    const double beta[4] = { 16.45934362, 299.1326597, 7833.604578, 420122.2104 };
    const double x3 = 222.0*(beta[0] - (1/222.0)*x[0] + (1/222.0)*x[1] + (1/222.0)*x[2])/219.0;
	double y2 = (1/222.0)*x[0]*x[0] + (1/222.0)*x[1]*x[1] + (1/222.0)*x[2]*x[2] + (219/222.0)*x3*x3 - beta[1];
	double y3 = (1/222.0)*x[0]*x[0]*x[0] + (1/222.0)*x[1]*x[1]*x[1] + (1/222)*x[2]*x[2]*x[2] + (219/222.0)*x3*x3*x3 - beta[2];
	double y4 = (1/222.0)*x[0]*x[0]*x[0]*x[0] + (1/222.0)*x[1]*x[1]*x[1]*x[1] + (1/222.0)*x[2]*x[2]*x[2]*x[2] + (219/222.0)*x3*x3*x3*x3 - beta[3];

    return sqrt(y2*y2+y3*y3+y4*y4);
}

int mainx(int argc, char* argv[]) {
	nlopt_opt opt;
    
    //finished at f(78.0242,16.5017,-7.84708) = 0.0146327154
    //finished at f(16.306,78.028,-0.697541) = 0.1118539931
    //finished at f(78.0261,16.46,5.49311) = 3.887896043e-05

    const double lb[3] = { 0, 0, 0 };
    
	opt = nlopt_create(NLOPT_LN_COBYLA, 3); /* algorithm and dimensionality */
	nlopt_set_min_objective(opt, hfunc, NULL);
	nlopt_set_xtol_rel(opt, 1e-7);
    nlopt_set_lower_bounds(opt, lb);
    
    // finished at f(93.4599,26.1291,23.8877) = 118886.9065


	double x[3] = { 90, 34.5, 20.625 };  /* some initial guess */
    cout << "initial f: " << hfunc(3, x, NULL, NULL) << endl;
	double minf; /* the minimum objective value, upon return */
    nlopt_result res = nlopt_optimize(opt, x, &minf);
	if (res < 0) {
		printf("nlopt failed! (%d)\n", res);
    }
    printf("finished at f(%g,%g,%g) = %0.10g\n", x[0], x[1], x[2],minf);
	
	nlopt_destroy(opt);
    return 0;
}

int main(int argc, char* argv[]) {
    for (double j = 34; j < 36; j += .01) {
        double x[3] = {93.15, j, 20.625 };
        double y = hfunc(3, x, NULL, NULL);
        cout << j << "\t" << y << endl;
    }
    return 0;
}