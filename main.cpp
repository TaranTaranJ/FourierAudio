#include <SDL2/SDL.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>

// Application window dimensions
const int WINDOW_WIDTH = 1800;
const int WINDOW_HEIGHT = 600;

using complexArr = std::complex<double>;

class PCMVisualizer {
public:
    PCMVisualizer(const char* filename);
    ~PCMVisualizer();
    void run();

private:
    void loadPCMData(const char* filename);
    void renderWaveform();
    void renderFrequencyDomain();
    void computeFFT();
    void FFT1D(std::vector<complexArr>& x);

    SDL_Window* m_window;
    SDL_Renderer* m_renderer;
    std::vector<short> m_pcmSamples;
    std::vector<complexArr> m_fftResult;
    bool displayFrequencyDomain;
};

PCMVisualizer::PCMVisualizer(const char* filename) : displayFrequencyDomain(false) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "Failed to initialize SDL: " << SDL_GetError() << std::endl;
        exit(1);
    }

    m_window = SDL_CreateWindow("PCM Waveform and Frequency Domain Visualization",
                                SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_SHOWN);

    if (!m_window) {
        std::cerr << "Failed to create window: " << SDL_GetError() << std::endl;
        SDL_Quit();
        exit(1);
    }

    m_renderer = SDL_CreateRenderer(m_window, -1, SDL_RENDERER_ACCELERATED);
    if (!m_renderer) {
        std::cerr << "Failed to create renderer: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(m_window);
        SDL_Quit();
        exit(1);
    }

    loadPCMData(filename);
    computeFFT();
}

PCMVisualizer::~PCMVisualizer() {
    SDL_DestroyRenderer(m_renderer);
    SDL_DestroyWindow(m_window);
    SDL_Quit();
}

void PCMVisualizer::loadPCMData(const char* filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    std::vector<char> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    m_pcmSamples.resize(buffer.size() / 2);
    for (size_t i = 0; i < buffer.size(); i += 2) {
        m_pcmSamples[i / 2] = static_cast<short>(buffer[i] | (buffer[i + 1] << 8));
    }

    if (m_pcmSamples.empty()) {
        std::cerr << "PCM data is empty!" << std::endl;
        exit(1);
    }

    std::cout << "Loaded " << m_pcmSamples.size() << " PCM samples." << std::endl;
}

void PCMVisualizer::FFT1D(std::vector<complexArr>& x) {
    const size_t N = x.size();
    if (N <= 1) return;

    std::vector<complexArr> even(N / 2);
    std::vector<complexArr> odd(N / 2);
    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    FFT1D(even);
    FFT1D(odd);

    for (size_t k = 0; k < N / 2; ++k) {
        complexArr t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}

void PCMVisualizer::computeFFT() {
    size_t N = m_pcmSamples.size();
    std::vector<complexArr> complexSamples(N);

    for (size_t i = 0; i < N; ++i) {
        complexSamples[i] = static_cast<double>(m_pcmSamples[i]);
    }

    FFT1D(complexSamples);
    m_fftResult = complexSamples;
    std::cout << "Computed FFT of " << N << " samples." << std::endl;
}

void PCMVisualizer::renderWaveform() {
    SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
    SDL_RenderClear(m_renderer);

    SDL_SetRenderDrawColor(m_renderer, 0, 255, 0, 255);

    int mid = WINDOW_HEIGHT / 2;
    size_t step = m_pcmSamples.size() / WINDOW_WIDTH;

    int x1 = 0;
    int sampleIndex = x1 * step;
    int sample = m_pcmSamples[sampleIndex] / (32768 / mid);
    int y1 = mid - sample;

    for (int x2 = 1; x2 < WINDOW_WIDTH; ++x2) {
        sampleIndex = x2 * step;
        sample = m_pcmSamples[sampleIndex] / (32768 / mid);
        int y2 = mid - sample;

        SDL_RenderDrawLine(m_renderer, x1, y1, x2, y2);

        x1 = x2;
        y1 = y2;
    }

    SDL_RenderPresent(m_renderer);
}

void PCMVisualizer::renderFrequencyDomain() {
    SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
    SDL_RenderClear(m_renderer);

    SDL_SetRenderDrawColor(m_renderer, 0, 0, 255, 255);

    int mid = WINDOW_HEIGHT / 2;
    size_t step = m_fftResult.size() / WINDOW_WIDTH;

    int x1 = 0;
    int sampleIndex = x1 * step;
    double magnitudeValue = std::abs(m_fftResult[sampleIndex]) / (m_fftResult.size() / mid);
    int y1 = mid - static_cast<int>(magnitudeValue);

    for (int x2 = 1; x2 < WINDOW_WIDTH; ++x2) {
        sampleIndex = x2 * step;
        magnitudeValue = std::abs(m_fftResult[sampleIndex]) / (m_fftResult.size() / mid);
        int y2 = mid - static_cast<int>(magnitudeValue);

        SDL_RenderDrawLine(m_renderer, x1, y1, x2, y2);

        x1 = x2;
        y1 = y2;
    }

    SDL_RenderPresent(m_renderer);
}

void PCMVisualizer::run() {
    bool running = true;
    SDL_Event event;

    Uint32 startTime = SDL_GetTicks();

    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }

        Uint32 currentTime = SDL_GetTicks();
        if (currentTime - startTime < 3000) {
            renderWaveform();
        } else {
            renderFrequencyDomain();
        }

        SDL_Delay(16); // Approx. 60 FPS
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pcm file>" << std::endl;
        return 1;
    }

    PCMVisualizer visualizer(argv[1]);
    visualizer.run();

    return 0;
}
