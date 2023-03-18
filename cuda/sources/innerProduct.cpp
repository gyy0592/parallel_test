#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <iostream>
#include <vector>
#include <chrono>

#define N 4000 // Size of vectors

// Kernel function to compute the inner product of two vectors on GPU
__global__ void innerProductKernel(float* A, float* B, float* result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        result[idx] = A[idx] * B[idx];
    }
}

// Function to compute the inner product of two vectors using CUDA
float innerProductCUDA(float* A, float* B) {
    float* d_A;
    float* d_B;
    float* d_result;
    float result = 0;

    cudaMalloc(&d_A, N * sizeof(float));
    cudaMalloc(&d_B, N * sizeof(float));
    cudaMalloc(&d_result, N * sizeof(float));

    cudaMemcpy(d_A, A, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, N * sizeof(float), cudaMemcpyHostToDevice);

    dim3 dimBlock(256, 1, 1);
    dim3 dimGrid((N + dimBlock.x - 1) / dimBlock.x, 1, 1);

    innerProductKernel<<<dimGrid, dimBlock>>>(d_A, d_B, d_result);

    cudaMemcpy(&result, d_result, sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_result);

    return result;
}

// Function to compute the inner product of two vectors using the CPU
float innerProductCPU(float* A, float* B) {
    float result = 0;
    for (int i = 0; i < N; i++) {
        result += A[i] * B[i];
    }
    return result;
}

int main() {
    std::vector<float> A(N);
    std::vector<float> B(N);

    // Initialize vectors with random values
    for (int i = 0; i < N; i++) {
        A[i] = rand() % 10;
        B[i] = rand() % 10;
    }

    // Compute inner product using CUDA
    auto start_cuda = std::chrono::high_resolution_clock::now();
    float result_cuda = innerProductCUDA(A.data(), B.data());
    auto stop_cuda = std::chrono::high_resolution_clock::now();
    auto duration_cuda = std::chrono::duration_cast<std::chrono::microseconds>(stop_cuda - start_cuda);

    // Compute inner product using CPU
    auto start_cpu = std::chrono::high_resolution_clock::now();
    float result_cpu = innerProductCPU(A.data(), B.data());
    auto stop_cpu = std::chrono::high_resolution_clock::now();
    auto duration_cpu = std::chrono::duration_cast<std::chrono::microseconds>(stop_cpu - start_cpu);

    // Check if results are equal
    if (result_cpu == result_cuda) {
        std::cout << "Results are equal\n";
    } else {
        std::cout << "Results are not equal\n";
    }

    // Print execution times
    std::cout << "CUDA: " << duration_cuda.count() << " microseconds\n";
    std::cout << "CPU: " << duration_cpu.count() << " microseconds\n";

    return 0;
}
