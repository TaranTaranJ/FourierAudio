#include <SDL2/SDL.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>

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
    void computeFFT();
    void FFT1D(std::vector<complexArr>& x);
    void renderWaveform();
    void renderFrequencyDomain();

    SDL_Window* m_window;
    SDL_Renderer* m_renderer;
    std::vector<short> m_pcmSamples;
    std::vector<complexArr> m_fftResult;
};

PCMVisualizer::PCMVisualizer(const char* filename) {
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

    // Determine the file size
    file.seekg(0, std::ios::end);
    std::streampos fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    // Read the entire file into a buffer
    std::vector<char> buffer(fileSize);
    if (!file.read(buffer.data(), fileSize)) {
        std::cerr << "Error reading file: " << filename << std::endl;
        file.close();
        exit(1);
    }
    file.close();

    // Convert the buffer to PCM samples
    size_t numSamples = fileSize / 2; // Each sample is 2 bytes (16 bits)
    m_pcmSamples.resize(numSamples);

    for (size_t i = 0; i < numSamples; ++i) {
        m_pcmSamples[i] = static_cast<short>((buffer[i * 2] & 0xFF) | (buffer[i * 2 + 1] << 8));
    }

    if (m_pcmSamples.empty()) {
        std::cerr << "PCM data is empty!" << std::endl;
        exit(1);
    }

    std::cout << "Loaded " << m_pcmSamples.size() << " PCM samples." << std::endl;
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

void PCMVisualizer::renderWaveform() {
    SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
    SDL_RenderClear(m_renderer);

    SDL_SetRenderDrawColor(m_renderer, 0, 255, 0, 255);

    int mid = WINDOW_HEIGHT / 2;
    size_t numSamples = m_pcmSamples.size();
    float step = static_cast<float>(numSamples) / static_cast<float>(WINDOW_WIDTH);

    int x1 = 0;
    int sampleIndex = 0;
    float scale = static_cast<float>(mid) / 32768.0f; // Scale factor for y-axis

    int y1 = mid - static_cast<int>(m_pcmSamples[sampleIndex] * scale);

    for (int x2 = 1; x2 < WINDOW_WIDTH; ++x2) {
        sampleIndex = static_cast<int>(x2 * step);
        if (sampleIndex >= numSamples) break; // Ensure we don't access out of bounds

        int y2 = mid - static_cast<int>(m_pcmSamples[sampleIndex] * scale);

        SDL_RenderDrawLine(m_renderer, x1, y1, x2, y2);

        x1 = x2;
        y1 = y2;
    }

    SDL_RenderPresent(m_renderer);
}




void PCMVisualizer::renderFrequencyDomain() {
    SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
    SDL_RenderClear(m_renderer);

    SDL_SetRenderDrawColor(m_renderer, 0, 255, 0, 255);

    int mid = WINDOW_HEIGHT / 2;
    size_t numBins = m_fftResult.size() / 2;
    float x_scale = static_cast<float>(WINDOW_WIDTH) / static_cast<float>(numBins);

    float maxMagnitude = 0.0f;
    for (size_t i = 0; i < numBins; ++i) {
        float magnitude = std::abs(m_fftResult[i]);
        if (magnitude > maxMagnitude) {
            maxMagnitude = magnitude;
        }
    }

    float y_scale = static_cast<float>(mid) / maxMagnitude;

    for (size_t i = 1; i < numBins; ++i) {
        int x1 = static_cast<int>((i - 1) * x_scale);
        int y1 = mid - static_cast<int>(std::abs(m_fftResult[i - 1]) * y_scale);

        int x2 = static_cast<int>(i * x_scale);
        int y2 = mid - static_cast<int>(std::abs(m_fftResult[i]) * y_scale);

        // Ensure y1 and y2 are within bounds
        y1 = std::max(0, std::min(WINDOW_HEIGHT, y1));
        y2 = std::max(0, std::min(WINDOW_HEIGHT, y2));

        SDL_RenderDrawLine(m_renderer, x1, y1, x2, y2);
    }

    SDL_RenderPresent(m_renderer);
}




void PCMVisualizer::run() {
    bool running = true;
    SDL_Event event;

    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }

        renderWaveform();
        SDL_Delay(3000);
        renderFrequencyDomain();
        SDL_Delay(3000); // Approx. 60 FPS
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
