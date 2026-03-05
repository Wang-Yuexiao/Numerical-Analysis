#include <stdio.h>
#include <math.h>

float get_eps_float() {
    float eps = 1.0f;
    while (1.0f + eps / 2.0f > 1.0f) {
        eps /= 2.0f;
    }
    return eps;
}

double get_eps_double() {
    double eps = 1.0;
    while (1.0 + eps / 2.0 > 1.0) {
        eps /= 2.0;
    }
    return eps;
}

int main() {
    float float_eps = get_eps_float();
    double double_eps = get_eps_double();

    printf("Computed float epsilon: %.8e\n", float_eps);
    printf("Computed double epsilon: %.16e\n", double_eps);

    return 0;
}
