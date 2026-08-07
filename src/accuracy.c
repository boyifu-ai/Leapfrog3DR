#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

static FILE *open_file_or_die(const char *path);
static int count_values_in_file(FILE *file);
static void print_usage(const char *prog);

int main(int argc, char **argv)
{
    FILE *fp_r;
    FILE *fp_phi2;
    FILE *fp_phi4;
    FILE *fp_phi8;
    FILE *fp_phi16;
    FILE *fp_phi32;
    int lines_r;
    int lines_phi2;
    int lines_phi4;
    int lines_phi8;
    int lines_phi16;
    int lines_phi32;
    int i;
    double err_2 = 0.0;
    double err_4 = 0.0;
    double err_8 = 0.0;
    double err_16 = 0.0;
    double err_32 = 0.0;
    double value_r;
    double value_phi;
    double h = 1.0 / (double)gnx;

    if (argc != 7)
    {
        print_usage(argv[0]);
        return 2;
    }

    fp_r = open_file_or_die(argv[1]);
    fp_phi2 = open_file_or_die(argv[2]);
    fp_phi4 = open_file_or_die(argv[3]);
    fp_phi8 = open_file_or_die(argv[4]);
    fp_phi16 = open_file_or_die(argv[5]);
    fp_phi32 = open_file_or_die(argv[6]);

    lines_r = count_values_in_file(fp_r);
    lines_phi2 = count_values_in_file(fp_phi2);
    lines_phi4 = count_values_in_file(fp_phi4);
    lines_phi8 = count_values_in_file(fp_phi8);
    lines_phi16 = count_values_in_file(fp_phi16);
    lines_phi32 = count_values_in_file(fp_phi32);

    printf("File reference has %d values.\n", lines_r);
    printf("File phi2 has %d values.\n", lines_phi2);
    printf("File phi4 has %d values.\n", lines_phi4);
    printf("File phi8 has %d values.\n", lines_phi8);
    printf("File phi16 has %d values.\n", lines_phi16);
    printf("File phi32 has %d values.\n", lines_phi32);

    if (lines_r != lines_phi2 || lines_r != lines_phi4 || lines_r != lines_phi8 ||
        lines_r != lines_phi16 || lines_r != lines_phi32)
    {
        fprintf(stderr, "Error: all input files must contain the same number of values.\n");
        return 1;
    }

    for (i = 0; i < lines_r; i++)
    {
        fscanf(fp_r, "%lf", &value_r);
        fscanf(fp_phi2, "%lf", &value_phi);
        err_2 += pow(value_phi - value_r, 2);
    }

    rewind(fp_r);
    for (i = 0; i < lines_r; i++)
    {
        fscanf(fp_r, "%lf", &value_r);
        fscanf(fp_phi4, "%lf", &value_phi);
        err_4 += pow(value_phi - value_r, 2);
    }

    rewind(fp_r);
    for (i = 0; i < lines_r; i++)
    {
        fscanf(fp_r, "%lf", &value_r);
        fscanf(fp_phi8, "%lf", &value_phi);
        err_8 += pow(value_phi - value_r, 2);
    }

    rewind(fp_r);
    for (i = 0; i < lines_r; i++)
    {
        fscanf(fp_r, "%lf", &value_r);
        fscanf(fp_phi16, "%lf", &value_phi);
        err_16 += pow(value_phi - value_r, 2);
    }

    rewind(fp_r);
    for (i = 0; i < lines_r; i++)
    {
        fscanf(fp_r, "%lf", &value_r);
        fscanf(fp_phi32, "%lf", &value_phi);
        err_32 += pow(value_phi - value_r, 2);
    }

    err_2 = sqrt(h * h * h * err_2);
    err_4 = sqrt(h * h * h * err_4);
    err_8 = sqrt(h * h * h * err_8);
    err_16 = sqrt(h * h * h * err_16);
    err_32 = sqrt(h * h * h * err_32);

    printf("h: %f\n", h);
    printf("Error_2: %e\n", err_2);
    printf("Error_4: %e\n", err_4);
    printf("Error_8: %e\n", err_8);
    printf("Error_16: %e\n", err_16);
    printf("Error_32: %e\n", err_32);
    printf("log_1: %f\n", log(err_4 / err_2) / log(2.0));
    printf("log_2: %f\n", log(err_8 / err_4) / log(2.0));
    printf("log_3: %f\n", log(err_16 / err_8) / log(2.0));
    printf("log_4: %f\n", log(err_32 / err_16) / log(2.0));

    fclose(fp_r);
    fclose(fp_phi2);
    fclose(fp_phi4);
    fclose(fp_phi8);
    fclose(fp_phi16);
    fclose(fp_phi32);
    return 0;
}

static FILE *open_file_or_die(const char *path)
{
    FILE *file = fopen(path, "r");
    if (file == NULL)
    {
        fprintf(stderr, "Error: could not open %s\n", path);
        exit(1);
    }
    return file;
}

static int count_values_in_file(FILE *file)
{
    int count = 0;
    double value;

    while (fscanf(file, "%lf", &value) == 1)
    {
        count++;
    }
    rewind(file);
    return count;
}

static void print_usage(const char *prog)
{
    fprintf(stderr, "Usage: %s reference phi_dt2 phi_dt4 phi_dt8 phi_dt16 phi_dt32\n", prog);
}
