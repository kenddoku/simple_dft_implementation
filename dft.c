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

// Since I dont use any analog devices I simulate analog signal by generating
// much more samples of desired curve (combination of sine waves)
int generate_signal(
    double* signal, size_t signal_len, double sampling_time,
    double sine_freqs[], double initial_phases[], size_t num_freqs
) {
    for(size_t freq_index=0; freq_index<num_freqs; freq_index++)
        for(size_t signal_point=0; signal_point<signal_len; signal_point++) {
            double t = signal_point * sampling_time / signal_len;
            signal[signal_point] += sin(2*M_PI * sine_freqs[freq_index] * t + initial_phases[freq_index]);
    }

    FILE *fp = fopen("signal.txt", "w");
    if(fp == NULL) {
        perror("fopen");
        return -1;
    }
    fprintf(fp, "sine_freqs[Hz]: ");
    for(size_t freq_index=0; freq_index<num_freqs; freq_index++) {
        fprintf(fp, "%.6f ", sine_freqs[freq_index]);
    }
    fprintf(fp, "\ninitial_phases[rad]: ");
    for(size_t freq_index=0; freq_index<num_freqs; freq_index++) {
        fprintf(fp, "%.6f ", initial_phases[freq_index]);
    }
    fprintf(fp, "\n");
    for(size_t signal_point=0; signal_point<signal_len; signal_point++) {
        double t = signal_point * sampling_time / signal_len;
        fprintf(fp, "%.6f %.6f\n", t, signal[signal_point]);
    }

    fclose(fp);

    return 0;
}

int generate_sampled_signal(
    double* samples, double sampling_freq, double samples_len,
    double sine_freqs[], double initial_phases[], size_t num_freqs
) {
    for(size_t i=0; i<samples_len; i++) {
        double t = i / sampling_freq;
        for(size_t freq_index=0; freq_index<num_freqs; freq_index++) {
            samples[i] += sin(2*M_PI * sine_freqs[freq_index] * t + initial_phases[freq_index]);
        }
    }

    FILE *fp = fopen("samples.txt", "w");
    if(fp == NULL) {
        perror("fopen");
        return -1;
    }
    
    for(size_t i=0; i<samples_len; i++) {
        double t = i / sampling_freq;
        fprintf(fp, "%.6f %.6f\n", t, samples[i]);
    }

    return 0;
    // DOKOŃCZ ZAPIS DO PLIKU I ZOBACZ CZY DZIALA BO NIE MA SENSU MIEC SYGNALU JAK NIE JEST ANALOGOWY-
}


int sample_signal(
    double *signal, size_t signal_len,
    double *samples, size_t samples_len,
    double sampling_time
) {
    size_t signal_sample_len_ratio = signal_len / samples_len; 
    for(size_t i=0; i<samples_len; i++) {
        size_t sampling_point = i * (signal_sample_len_ratio);
        samples[i] = signal[sampling_point];
    }

    FILE *fp = fopen("samples.txt", "w");
    if(fp == NULL) {
        perror("fopen");
        return -1;
    }
    fprintf(fp, "sampling_time[s]: %f\n", sampling_time);
    for(size_t i=0; i<samples_len; i++) {
        fprintf(fp, "%.6f\n", samples[i]);
    }

    fclose(fp);

    return 0;
}

int main() {

    double SAMPLING_FREQ = 1024.0;   // [Hz]
    double SAMPLING_TIME = 1.0;    // [s]
    size_t SAMPLES_LEN = SAMPLING_FREQ * SAMPLING_TIME;
    size_t SIGNAL_LEN  = 10 * SAMPLES_LEN;
    double NYQUIST_FREQ = SAMPLES_LEN / 2.0; // [Hz] max frequency that can be reliably extracted from the signal
    
    // Specifying frequency components of a signal ================================================
    double sine_freqs[]     = {2.3, 5.7, 11.2, 27.77}; // [Hz]  
    double initial_phases[] = {M_PI/7.0, M_PI/2.3,  0.0,  11.7}; // [rad]
    size_t num_freqs        = sizeof(sine_freqs)/sizeof(sine_freqs[0]);

    // Allocating memory for signal ===============================================================
    double *signal = (double*)calloc(SIGNAL_LEN, sizeof(double));
    // Allocating memory for samples vector and zero initializing it ==============================
    double *samples = (double*)calloc(SAMPLES_LEN, sizeof(double));

    // Generating desired signal ==================================================================
    generate_signal(signal, SIGNAL_LEN, SAMPLING_TIME, sine_freqs, initial_phases, num_freqs);
    // Sampling the signal
    // sample_signal(signal, SIGNAL_LEN, samples, SAMPLES_LEN, SAMPLING_TIME);
    generate_sampled_signal(samples, SAMPLING_FREQ, SAMPLES_LEN, sine_freqs, initial_phases, num_freqs);

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