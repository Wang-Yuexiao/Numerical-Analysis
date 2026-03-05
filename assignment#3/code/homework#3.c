#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "nrutil.h"

#define NRANSI
#define PRINT_EXTRA 0
#define NP 20
#define MP 20
#define MAXSTR 80

typedef float const*const*const NonModifiableMatrix;
typedef float const*const NonModifiableVector;

void mprove(float **a, float **alud, int n, int indx[], float b[], float x[]);
void lubksb(float **a, int n, int *indx, float b[]);
void gaussj(float **a, int n, float **b, int m);
void ludcmp(float **a, int n, int *indx, float *d);
void svdcmp(float **a, int m, int n, float w[], float **v);

int read_matrix(FILE* fp, int* n, int* m, float** a, float* b);
void hw1_gaussj(int m, int n, NonModifiableMatrix a, NonModifiableVector b);
void hw1_LU(int m, int n, NonModifiableMatrix a, NonModifiableVector b);
void hw1_SVD(int m, int n, NonModifiableMatrix a, NonModifiableVector b);
void hw2_3(int m, int n, NonModifiableMatrix a, NonModifiableVector b);

int main(int argc, char const *argv[])
{
    int n,m,err;
	float **a,*b;
    const char *filename;
	
	FILE *fp;

    if(argc != 2){
        printf("Usage: %s <input file>\n", argv[0]);
        return -1;
    }
    filename = argv[1];

    a=matrix(1,MP,1,NP);
    b=vector(1,NP);
	
    if ((fp = fopen(filename, "r")) == NULL){
        char buf[MAXSTR];
        sprintf(buf, "Data file %s not found\n", filename);
		nrerror(buf);
        exit(-1);
    }

    if((err = read_matrix(fp, &m, &n, a, b))){
        char buf[MAXSTR];
        sprintf(buf, "Error reading matrix a: %d\n", err);
        nrerror(buf);
        exit(-1);
    }

    //lineq1 can't use Guass-Jordan Elimination
    hw1_gaussj(m, n, a, b);
    hw1_LU(m, n, a, b);
    hw1_SVD(m, n, a, b);
    hw2_3(m, n, a, b);

    fclose(fp);
    
	free_matrix(a,1,MP,1,NP);
    free_vector(b,1,NP);

    return 0;
}

int read_matrix(FILE* fp, int* m, int* n, float** a, float* b){
    int k,l;

    if(feof(fp))
        return -1;
	
    if(fscanf(fp,"%d %d ", m, n) != 2)
        return -2;

    for (k=1;k<=*m;k++){
        for (l=1;l<=*n;l++) {
            if(fscanf(fp,"%f ", &a[k][l]) != 1)
                return -3;
        }
    }

    for (k=1;k<=*n;k++){
        if(fscanf(fp,"%f ", &b[k]) != 1)
            return -3;
    }

    return 0;
}

void hw1_gaussj(int m, int n, NonModifiableMatrix a, NonModifiableVector b){
    int l,k,i,j;
    float **a_inv,**u,**x;

    assert(m == n);

    a_inv=matrix(1,NP,1,NP);
	u=matrix(1,NP,1,NP);
    x=matrix(1,NP,1,NP);

    printf("Using Gauss-Jordan Elimination\n");

    /* save matrices for later testing of results */
    for (l=1;l<=n;l++) {
        for (k=1;k<=n;k++){
            a_inv[k][l]=a[k][l];
        }
    }

    /* invert matrix */
    gaussj(a_inv,n,x,n);

    printf("Solution vector for the equations:\n");
    for (k=1;k<=n;k++) {
        float v = 0.0;
        for(l=1;l<=n;l++){
            v += a_inv[k][l]*b[l];
        }
        printf("%12.6f", v);
    }
    printf("\n");
    printf("------------------------------------------------------------------\n");

    free_matrix(x,1,NP,1,NP);
	free_matrix(u,1,NP,1,NP);
	free_matrix(a_inv,1,NP,1,NP);
}

void hw1_LU(int m, int n, NonModifiableMatrix a, NonModifiableVector b){
    int l,k,i,j,dum;
    float d;
    int *indx, *jndx;
    float **aa, **xl,**xu,**x,*ans;

    indx=ivector(1,NP);
	jndx=ivector(1,NP);

    aa=matrix(1,MP,1,NP);
    xl=matrix(1,NP,1,NP);
	xu=matrix(1,NP,1,NP);
	x=matrix(1,NP,1,NP);
    ans=vector(1,MP);

    printf("Using LU Decomposition\n");

    for (k=1;k<=m;k++) {
        for(l=1;l<=n;l++){
            aa[k][l]=a[k][l];
        }
        ans[k]=b[k];
    }

    ludcmp(aa,n,indx,&d);
    
    lubksb(aa,n,indx,ans);

    printf("Solution vector for the equations:\n");
    for (k=1;k<=m;k++) printf("%12.6f",ans[k]);
    printf("\n");

    printf("------------------------------------------------------------------\n");

    free_vector(ans,1,MP);

    free_ivector(indx, 1, NP);
	free_ivector(jndx, 1, NP);

    free_matrix(aa, 1,MP,1,NP);
    free_matrix(xl, 1, NP, 1, NP);
	free_matrix(xu, 1, NP, 1, NP);
	free_matrix(x, 1, NP, 1, NP);
}

void hw1_SVD(int m, int n, NonModifiableMatrix a, NonModifiableVector b){
    int j,k,l;
	float **aa,*w,**u,**v,*bb;

    aa=matrix(1,NP,1,MP);
    w=vector(1,NP);
    bb=vector(1,MP);
	u=matrix(1,NP,1,MP);
	v=matrix(1,NP,1,NP);

    printf("Using Singular Value Decomposition\n");

    /* copy original matrix into u */
    for (k=1;k<=m;k++){
        bb[k] = b[k];
        for (l=1;l<=n;l++) {
            aa[k][l]=a[k][l];
            u[k][l]=a[k][l];
        }
    }

    /* perform decomposition */
    svdcmp(u,m,m,w,v);

    printf("Solution vector for the equations:\n");
    for (k=1;k<=n;k++) {
        const float w_inv = (fabs(w[k]) >= 0.5e-6 ? 1.0 / w[k] : 0.0);
        bb[k] = 0.0;
        for(j=1;j<=m;j++)
            bb[k] += u[j][k] * b[j];
        bb[k] *= w_inv;
    }
    for (k=1;k<=n;++k){
        float ans = 0.0;
        for (j=1;j<=n;++j)
            ans += v[k][j] * bb[j];
        printf("%12.6f", ans);
    }

    printf("\n");

    printf("------------------------------------------------------------------\n");

    free_matrix(aa, 1, NP, 1, MP);
    free_vector(w, 1, NP);
    free_vector(bb, 1, MP);
    free_matrix(u, 1, NP, 1, MP);
    free_matrix(v, 1, NP, 1, NP);
}

void hw2_3(int m, int n, NonModifiableMatrix a, NonModifiableVector b){
    int i,j,*indx;
	float d,*x,**aa,*col,**a_inv;

    printf("Using the iterative improvement mprove method\n");

	indx=ivector(1,n);
	x=vector(1,n);
	aa=matrix(1,n,1,n);
    a_inv=matrix(1,n,1,n);
    col=vector(1,n);
	for (i=1;i<=n;i++) {
		x[i]=b[i];
		for (j=1;j<=n;j++)
			aa[i][j]=a[i][j];
	}
	ludcmp(aa,n,indx,&d);
	lubksb(aa,n,indx,x);

	printf("Solution vector for the equations:\n");
	for (i=1;i<=n;i++) printf("%12.6f",x[i]);
	printf("\n");

    for (i=1;i<=n;i++) {
        for (j=1;j<=n;j++) col[j] = 0.0;
        col[i] = 1.0;
        lubksb(aa,n,indx,col);
        for (j=1;j<=n;j++) a_inv[j][i] = col[j];
    }
    printf("Inverse of the matrix:\n");
    for (i=1;i<=n;i++) {
        for (j=1;j<=n;j++) printf("%12.6f",a_inv[i][j]);
        printf("\n");
    }

    for (i=1;i<=m;i++)
        d *= aa[i][i];
    printf("Determinant of the matrix: %12.6f\n", d);

	free_matrix(aa,1,n,1,n);
	free_vector(x,1,n);
	free_ivector(indx,1,n);
    free_matrix(a_inv,1,n,1,n);
    free_vector(col,1,n);
}

#undef NRANSI
