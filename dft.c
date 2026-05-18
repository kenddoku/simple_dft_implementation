// This program serves as a simple implementation of the Discrete Fourier Transform

// Discrete Fourier Transform can be calculated using expression below
// f_j_hat = sum_k=0^(n-1){f_k * exp{-2pi * i * j/n * k}}
// where f_k is the kth data sample and f_k_hat is the kth element of a vector in
// frequency domain 
// we define omega=exp{-2pi * i / n} as the fundamental frequency (frequency resolution)
// and notice that exp{2pi * i * j/n * k} = omega_j^k = omega^(j*k)
// f_j_hat = sum_k=0^(n-1){f_k * omega_j^k}

// The inverse of the Discrete Fourier Transform 
// has almost identical form (change in the sign of exponent)
// f_k = (1/n) * sum_k=0^(n-1){f_k_hat * exp(+2pi * i * k/n)}

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

int perform_dft(double complex **transform_matrix, double *samples, double complex *outputs, size_t N) {
    if(N <= 0.0) {
        return -1;
    }

    double complex omega = cexp(2*M_PI * I / N);

    for(size_t k=0; k<N; k++) {
        for(size_t n=0; n<N; n++) {
            transform_matrix[k][n] = cpow(omega, k*n);
            outputs[k] += transform_matrix[k][n] * samples[n];
        }
    }

    return 0;
}

void get_freq_magnitudes_and_phases(double complex *dft_outputs, double *freq_magnitudes, double* phases, size_t N) {
    for(size_t i=0; i<N; i++) {
        double real = creal(dft_outputs[i]);
        double imag = cimag(dft_outputs[i]);
        freq_magnitudes[i] = sqrt(real*real + imag*imag);
        phases[i] = carg(dft_outputs[i]);
    }
}

int save_dft_results(double sampling_freq, double* freq_magnitudes, double* phases, size_t N) {
    char filename[] = "freq_magnitudes_phases.txt";
    FILE *fp = fopen(filename, "w");

    if(fp == NULL) {
        perror("fopen");
        return -1;
    }

    // Saving only first half (up to Nyquist frequency)
    for(size_t i=0; i<= N/2; i++) {
        double freq = i * sampling_freq/N;
        fprintf(fp, "%.6f %.6f %.6f\n",
            freq,
            freq_magnitudes[i],
            phases[i]
        );
    }

    fclose(fp);

    return 0;
};

int main() {

    double SAMPLING_FREQ = 100.0;   // [Hz]
    double SAMPLING_TIME = 10;    // [s]
    size_t SAMPLES_LEN = SAMPLING_FREQ * SAMPLING_TIME;
    
    // Allocating memory for samples vector and zero initializing it ==============================
    double *samples = (double*)calloc(SAMPLES_LEN, sizeof(double));

    // Generating a simple sine wave (already sampled for simplicy; at least for now) =============
    // and saving it to a file
    double sine_freq = 10.0;
    double initial_phase = 0.0; // radians
    
    char samples_filename[] = "sampled_sine_wave.txt";
    FILE *samples_fp = fopen(samples_filename, "w");
    if(samples_fp == NULL) {
        perror("fopen");
        return -1;
    }

    fprintf(samples_fp, "sine_freq[Hz]: %.6f\n", sine_freq);
    fprintf(samples_fp, "initial_phase[rad]: %.6f\n", initial_phase);
    fprintf(samples_fp, "sampling_freq[Hz]: %.6f\n", SAMPLING_FREQ);
    fprintf(samples_fp, "sampling_time[s]: %.6f\n", SAMPLING_TIME);
    fprintf(samples_fp, "samples_len: %ld\n", SAMPLES_LEN);

    for(size_t i=0; i<SAMPLES_LEN; i++) {
        double freq_ratio = sine_freq / SAMPLING_FREQ;
        double t = i * SAMPLING_TIME/SAMPLES_LEN;
        samples[i] = sin(2*M_PI * i * freq_ratio + initial_phase);
        fprintf(samples_fp, "%.6f %.6f\n", t, samples[i]);
    }
    fclose(samples_fp);

    // Allocating memory for the Fourier Transform Matrix =========================================
    double complex **dft_matrix = (double complex**)malloc(SAMPLES_LEN * sizeof(double complex*));
    for(size_t i=0; i<SAMPLES_LEN; i++) {
        dft_matrix[i] = (double complex*)calloc(SAMPLES_LEN, sizeof(double complex));
    }
        
    // Allocating memory for DFT outputs ==========================================================
    double complex* outputs = (double complex*)calloc(SAMPLES_LEN, sizeof(double complex));

    perform_dft(dft_matrix, samples, outputs, SAMPLES_LEN);

    // Allocating memory for frequency magnitudes and phases ======================================
    double *freq_magnitudes = (double*)calloc(SAMPLES_LEN, sizeof(double));
    double *phases = (double*)calloc(SAMPLES_LEN, sizeof(double));

    get_freq_magnitudes_and_phases(outputs, freq_magnitudes, phases, SAMPLES_LEN);

    save_dft_results(SAMPLING_FREQ, freq_magnitudes, phases, SAMPLES_LEN);

    // Freeing memory =============================================================================
    free(phases);
    free(freq_magnitudes);
    free(outputs);
    free(samples);
    for(size_t i=0; i<SAMPLES_LEN; i++) {
        free(dft_matrix[i]);
    }
    free(dft_matrix);

    return 0;
}