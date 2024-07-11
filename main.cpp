#include <SDL2/SDL.h>
#include <iostream>
#include <fstream>
#include <vector>

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;

class PCMVisualizer {
public:
    PCMVisualizer(const char* filename);
    ~PCMVisualizer();
    void run();

private:
    void loadPCMData(const char* filename);
    void renderWaveform();

    SDL_Window* m_window;
    SDL_Renderer* m_renderer;
    std::vector<short> m_pcmSamples;
};

PCMVisualizer::PCMVisualizer(const char* filename) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "Failed to initialize SDL: " << SDL_GetError() << std::endl;
        exit(1);
    }

    m_window = SDL_CreateWindow("PCM Waveform Visualization",
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
}

void PCMVisualizer::renderWaveform() {
    SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
    SDL_RenderClear(m_renderer);

    SDL_SetRenderDrawColor(m_renderer, 0, 255, 0, 255);

    int mid = WINDOW_HEIGHT / 2;
    size_t step = m_pcmSamples.size() / WINDOW_WIDTH;

    // Initial point
    int x1 = 0;
    int sampleIndex = x1 * step;
    int sample = m_pcmSamples[sampleIndex] / (32768 / mid);
    int y1 = mid - sample;

    // Draw lines between successive points
    for (int x2 = 1; x2 < WINDOW_WIDTH; ++x2) {
        sampleIndex = x2 * step;
        sample = m_pcmSamples[sampleIndex] / (32768 / mid);
        int y2 = mid - sample;

        SDL_RenderDrawLine(m_renderer, x1, y1, x2, y2);

        // Update the previous point
        x1 = x2;
        y1 = y2;
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
