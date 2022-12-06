#include <stdio.h>

__global__ void dynamicReverse(int *d, int n) {
    extern __shared__ int s[];
    
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int t = threadIdx.x;
    int tr = 1024 - t - 1;

    s[t] = d[i];

    __syncthreads();

    d[i] = s[tr];
}

int main(void) {
    const int n = 10000;
    const int n_blocks = n / 1024 + (n % 1024 == 0 ? 0 : 1);

    printf("Numero di blocchi: %d\n", n_blocks);

    int a[n], r[n], d[n];

    for (int i = 0; i < n; i++) {
        a[i] = i;
        r[i] = n-i-1;
        d[i] = 0;
    }

    int *d_d;
    cudaMalloc(&d_d, n * sizeof(int)); 

    // run dynamic shared memory version
    cudaMemcpy(d_d, a, n*sizeof(int), cudaMemcpyHostToDevice);

    dynamicReverse<<<n_blocks, 1024, n*sizeof(int)>>>(d_d, n);

    cudaMemcpy(d, d_d, n * sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < n; i++) 
        if (d[i] != r[i]) printf("Error: d[%d]!=r[%d] (%d, %d)\n", i, i, d[i], r[i]);

    for (int i = 0; i < n; i++) 
        printf("d[%d] -> %d\n", i, d[i]);
}