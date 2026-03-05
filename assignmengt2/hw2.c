#include <math.h>
#include <stdio.h>
#include <sys/time.h>

#include "C:/c_project/NR_Assignment_2/NRs/ansi/recipes/bessj0.c"
#include "C:/c_project/NR_Assignment_2/NRs/ansi/recipes/bessj1.c"
#include "C:/c_project/NR_Assignment_2/NRs/ansi/recipes/zbrak.c"
#include "C:/c_project/NR_Assignment_2/NRs/ansi/recipes/rtbis.c"
#include "C:/c_project/NR_Assignment_2/NRs/ansi/recipes/rtflsp.c"
#include "C:/c_project/NR_Assignment_2/NRs/ansi/recipes/rtsec.c"
#include "C:/c_project/NR_Assignment_2/NRs/ansi/recipes/rtnewt.c"
#include "C:/c_project/NR_Assignment_2/NRs/ansi/recipes/rtsafe.c"
#include "C:/c_project/NR_Assignment_2/NRs/ansi/other/nrutil.c"

float myfunc(float);
void myfunc_(float, float *, float *);
void bessj0_(float, float *, float *);
void bessj0_(float, float *, float *);
float rtmuller(float (*func)(float), float x1, float x2, float xacc);

struct timeval start, end;

#define PRINT_AND_MEASURE(name, expr) \
	printf(name ": %f (Required time = %f)\n", (expr), num_steps);

int main(int argc, char const *argv[])
{
    static const float XACC = 1.0e-6;
    float xb1[101], xb2[101];
    int i,nb=20;
    double time_spent;


    float x1 = 1.0, x2 = 10.0;
    printf("Using bracketing method, find the possible roots of bessj0.\n");
    zbrak(bessj0, x1, x2, 100, xb1, xb2, &nb);

    //Bessel function
    //Bisection method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtbis(bessj0, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Bisection method root%d: %f\n",i, (rtbis(bessj0, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Linear Interpolation method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtflsp(bessj0, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Linear Interpolation root%d: %f\n",i, (rtflsp(bessj0, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Secant method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtsec(bessj0, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Secant method root%d: %f\n",i, (rtsec(bessj0, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Newton-Raphson method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtnewt(bessj0_, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Newton-Raphson method root%d: %f\n",i, (rtnewt(bessj0_, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Newton with bracketing
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtsafe(bessj0_, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Newton with bracketing root%d: %f\n",i, (rtsafe(bessj0_, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Muller Method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtmuller(bessj0, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Muller Method root%d: %f\n",i, (rtmuller(bessj0, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    printf("\n\n");
    //myfunction
	printf("Using bracketing method, find the possible roots of my nonlinear function f(x) = x^3 - 2x + 1\n");
	x1 = -5.0; x2 = 5.0;
    zbrak(myfunc, x1, x2, 100, xb1, xb2, &nb);
    //Bisection method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtbis(myfunc, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Bisection method root%d: %f\n",i, (rtbis(myfunc, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Linear Interpolation method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtflsp(myfunc, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Linear Interpolation root%d: %f\n",i, (rtflsp(myfunc, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Secant method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtsec(myfunc, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Secant method root%d: %f\n",i, (rtsec(myfunc, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Newton-Raphson method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtnewt(myfunc_, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Newton-Raphson method root%d: %f\n",i, (rtnewt(myfunc_, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Newton with bracketing
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtsafe(myfunc_, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Newton with bracketing root%d: %f\n",i, (rtsafe(myfunc_, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    //Muller Method
    gettimeofday(&start, NULL);
    for(int j = 0; j < 100000; j++){
        for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        rtmuller(myfunc, x1, x2, XACC);
        }
    }
    for (i=1;i<=nb;i++){
            x1 = xb1[i];
            x2 = xb2[i];
        printf("Muller Method root%d: %f\n",i, (rtmuller(myfunc, x1, x2, XACC)));
    }
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("(Require time = %f)\n", time_spent);

    return 0;
}

float myfunc(float x){
	// f(x) = x^3 - 2x + 1
	return x*x*x - 2*x + 1;
}

void myfunc_(float rtn, float *f, float *df){
	float x = rtn;
	*f = myfunc(x);
	// f'(x) = 3x^2 - 2
	*df = 3*x*x - 2;
}

float rtmuller(float (*func)(float), float x1, float x2, float xacc){
    float x3 = 0.5 * (x1 + x2);
    float x4;
    float a, b, c;
	const double eps = 1.084202e-19;
    do{
        c = func(x3);
		b = (x1 - x3) / (x2 - x3 + eps) / (x1 - x2 + eps) * (func(x2) - func(x3)) - (x2 - x3) / (x1 - x3 + eps) / (x1 - x2 + eps) * (func(x1) - func(x3));
        a = (func(x1) - func(x3)) / (x1 - x3 + eps) / (x1 - x2 + eps) - (func(x2) - func(x3)) / (x2 - x3 + eps) / (x1 - x2 + eps);
		if(b * b - 4 * a * c < 0){
			break;
		}
        x4 = x3 - 2 * c / (b + (b > 0 ? 1 : -1) * sqrt(b * b - 4 * a * c) + eps);
        x1 = x2; x2 = x3; x3 = x4;
    }while(fabs(func(x4)) > xacc);
    return x4;
}

void bessj0_(float rtn, float *f, float *df){
    *f = bessj0(rtn);
    *df = -bessj1(rtn);
}

