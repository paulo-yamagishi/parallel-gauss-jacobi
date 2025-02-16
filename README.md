# Parallel Gauss-Jacobi Solver

## Project Description

This project implements the **Gauss-Jacobi method** for solving systems of linear equations using **parallel programming** with POSIX threads (`pthread`). It compares the performance of the sequential and parallel implementations to analyze speedup and efficiency.

## Features

- **Sequential and Parallel Implementations**: Two versions of the Gauss-Jacobi algorithm for performance comparison.
- **Dynamic Matrix and Vector Generation**: Generates random matrices with diagonal dominance.
- **Multi-threaded Computation**: Utilizes `pthread` to distribute computation across multiple threads.
- **Performance Analysis**: Measures iteration count, execution time, speedup, and efficiency.

## Installation & Compilation

### Prerequisites

- **GCC Compiler**
- **POSIX Threads (pthread)**

### Compilation

To compile the program, run:

```sh
gcc -o gauss_jacobi_parallel parallel_Gauss_Jacobi.c -pthread -lm
```

## Usage

Run the program with:

```sh
./gauss_jacobi_parallel <matrix_size> <num_threads> <epsilon> <random_seed>
```

Where:

- `<matrix_size>`: Size of the square matrix (N × N)
- `<num_threads>`: Number of threads to use
- `<epsilon>`: Convergence threshold
- `<random_seed>`: Seed for random number generation

### Example

```sh
./gauss_jacobi_parallel 1000 4 0.0001 42
```

## Output

The program outputs:

- Number of iterations for convergence
- Execution time (sequential vs. parallel)
- Distance between sequential and parallel solutions
- Speedup and efficiency of the parallel implementation

## Authors

- **Paulo Yamagishi**
- **Gustavo Mariz** ([GitHub Profile](https://github.com/gustavomariz))

## License

This project is open-source and distributed under the **MIT License**.

