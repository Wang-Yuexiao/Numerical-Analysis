#include <stdio.h>
#include <float.h>

void machar(float *eps) {
    *eps = 1.0f;
    while (1.0f + *eps > 1.0f) {
        *eps /= 2.0f;
    }
    *eps *= 2.0f;
}

void machar_double(double *eps) {
    *eps = 1.0;
    while (1.0 + *eps > 1.0) {
        *eps /= 2.0;
    }
    *eps *= 2.0;
}

int main() {
    float float_eps;
    double double_eps;

    machar(&float_eps);
    machar_double(&double_eps);

    printf("float : %.8e\n", float_eps);
    printf("double : %.16e\n", double_eps);

    return 0;
}
