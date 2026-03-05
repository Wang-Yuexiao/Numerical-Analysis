#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define A -3
#define B 4
#define M 0.5
#define S 1.5

int sample_sizes[] = {100, 1000, 10000, 100000}; // 样本数量
int num_sample_sizes = sizeof(sample_sizes) / sizeof(sample_sizes[0]);

double generate_uniform(double a, double b) {
    return a + (b - a) * ((double)rand() / RAND_MAX);
}

double generate_gaussian(double mean, double stddev) {
    double u1 = ((double)rand() / RAND_MAX);
    double u2 = ((double)rand() / RAND_MAX);
    double z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
    return z0 * stddev + mean;
}

void generate_histogram(double (*rand_func)(double, double), double param1, double param2, int sample_size, const char *filename) {
    int bins[100] = {0}; // 用100个区间
    double min = -3;
    double max = 4;
    double interval = (max - min) / 100.0;

    // 生成样本并统计频率
    for (int i = 0; i < sample_size; i++) {
        double value = rand_func(param1, param2);
        
        if (value >= min && value < max) {
            int bin_index = (int)((value - min) / interval);
            bins[bin_index]++;
        }
    }

    // 将数据写入文件
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Failed to open file\n");
        return;
    }
    for (int i = 0; i < 100; i++) {
        fprintf(file, "%f %d\n", min + i * interval, bins[i]);
    }
    fclose(file);
}

int main() {
    srand(time(0)); // 初始化随机数种子

    for (int i = 0; i < num_sample_sizes; i++) {
        int sample_size = sample_sizes[i];

        // 生成均匀分布数据
        char uniform_filename[50];
        sprintf(uniform_filename, "uniform_%d.txt", sample_size);
        generate_histogram(generate_uniform, A, B, sample_size, uniform_filename);
        printf("Uniform distribution data saved to %s\n", uniform_filename);

        // 生成高斯分布数据
        char gaussian_filename[50];
        sprintf(gaussian_filename, "gaussian_%d.txt", sample_size);
        generate_histogram(generate_gaussian, M, S, sample_size, gaussian_filename);
        printf("Gaussian distribution data saved to %s\n", gaussian_filename);
    }

    return 0;
}