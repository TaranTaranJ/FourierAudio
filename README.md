#### Advanced Maths Resubmission 2024
# PCMVisualizer Code Breakdown
#
### How to Use:

1. Compilation: Ensure SDL2 library is linked (-lSDL2 for GCC).
2. Execution: ./PutPixel audio.pcm
   1. Replace audio.pcm with the files directory 
   2. Example: ./PutPixel /home/STUDENT ID/PutPixel/AUDIO-BEFORE.pcm
   
### Functionality:

This program visualizes PCM (Pulse Code Modulation)
audio data in both time domain and frequency domain using SDL
(Simple DirectMedia Layer) for graphics rendering.
It computes FFT (Fast Fourier Transform) to convert
PCM data from time domain to frequency domain and back.

1. Initialization: Initializes SDL for window creation and rendering.

```cpp
PCMVisualizer::PCMVisualizer(const char* filename)
        : m_targetFrequency(0.0), m_sampleRate(0.0), m_bandwidth(0.0) {
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
```
Once the code has been run and exit has been called.
```c++
PCMVisualizer::~PCMVisualizer() {
SDL_DestroyRenderer(m_renderer);
SDL_DestroyWindow(m_window);
SDL_Quit();
}
```
This closes and destroys all SDL processes.
#
2. PCM Loading: Loads PCM data from a binary file.
```cpp
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
    size_t numSamples = fileSize / sizeof(short); // Each sample is 2 bytes (16 bits)
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
```
This function encapsulates the process of reading PCM data from a binary file into memory.
It ensures proper error handling for file operations, calculates the necessary file size, 
reads the file contents into a buffer, converts the buffer contents into PCM samples,
and provides feedback on the number of samples loaded. 
This is crucial for subsequent processing and visualization of audio data within the PCMVisualizer application.
#
3. FFT Computation: Computes the FFT of loaded PCM data.
```cpp
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
```
This function initializes a vector of complex numbers from PCM samples, 
computes the one-dimensional Fast Fourier Transform (FFT) using FFT1D, and stores the result in m_fftResult. 
It prints a message indicating the number of samples processed.
#
```cpp
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
```
Implemented recursively, FFT1D performs the Cooley-Tukey algorithm to compute the FFT of a vector x containing complex numbers. 
It splits x into even and odd indexed elements, recursively computes FFT on these segments, applies the FFT butterfly operation to combine results, 
and transforms the input vector in-place.
#
```cpp
void PCMVisualizer::IFFT1D(std::vector<complexArr>& x) {
    // Conjugate the complex numbers
    for (auto& elem : x) {
        elem = std::conj(elem);
    }

    // Perform the forward FFT
    FFT1D(x);

    // Conjugate the complex numbers again
    for (auto& elem : x) {
        elem = std::conj(elem);
    }

    // Scale the numbers
    size_t N = x.size();
    for (auto& val : x) {
        val /= N;
    }

    std::cout << "Performed IFFT on " << N << " samples." << std::endl;
}
```
This function performs the Inverse Fast Fourier Transform (IFFT) on a vector x of complex numbers. It conjugates each element of x, 
performs FFT (which effectively computes IFFT due to conjugation), conjugates again, scales the results by dividing by the vector size, 
and prints a message indicating the number of samples processed.
#
4. Frequency Domain Manipulation: Removes a tone from the frequency domain.

```cpp
void PCMVisualizer::removeTone() {
size_t N = m_fftResult.size();
double binWidth = m_sampleRate / static_cast<double>(N); // Frequency resolution per bin

    // Calculate the index of the target frequency in the FFT result
    size_t targetIndex = static_cast<size_t>(m_targetFrequency / binWidth);
    size_t bandwidthBins = static_cast<size_t>(m_bandwidth / binWidth); // Number of bins to zero out around the target frequency

    // Ensure the targetIndex and bandwidthBins are within valid bounds
    if (targetIndex < N) {
        std::cout << "Target frequency: " << m_targetFrequency << " Hz" << std::endl;
        std::cout << "Bin width: " << binWidth << " Hz" << std::endl;
        std::cout << "Target index: " << targetIndex << std::endl;
        std::cout << "Bandwidth bins: " << bandwidthBins << std::endl;

        // Zero out the target frequency bin and its neighbors within the bandwidth
        for (size_t i = 0; i <= bandwidthBins; ++i) {
            if (targetIndex + i < N) {
                m_fftResult[targetIndex + i] = 0;
            }
            if (targetIndex >= i) {
                m_fftResult[targetIndex - i] = 0;
            }
        }

        // Since FFT result is symmetric, also zero out the corresponding negative frequency bins
        if (targetIndex != 0 && targetIndex != N - 1) {
            for (size_t i = 0; i <= bandwidthBins; ++i) {
                if (N - targetIndex + i < N) {
                    m_fftResult[N - targetIndex + i] = 0;
                }
                if (N - targetIndex >= i) {
                    m_fftResult[N - targetIndex - i] = 0;
                }
            }
        }
    } else {
        std::cerr << "Target index out of bounds: " << targetIndex << std::endl;
    }
}
```
The removeTone() function identifies the index (targetIndex) in the FFT result (m_fftResult) corresponding to a specified m_targetFrequency within a m_bandwidth. 
It then zeros out the magnitude values of the target frequency bin and its neighboring bins within the bandwidth.
This ensures that specific frequency components are effectively removed from the audio data in the frequency domain.
#

5. Conversion: Converts manipulated frequency domain back to time domain using IFFT.

```cpp
void PCMVisualizer::convertFrequencyToTimeDomain() {
IFFT1D(m_fftResult);

    // Convert the complex result back to PCM samples
    size_t N = m_fftResult.size();
    m_pcmSamples.resize(N);
    for (size_t i = 0; i < N; ++i) {
        m_pcmSamples[i] = static_cast<short>(std::real(m_fftResult[i]));
    }

    std::cout << "Converted frequency domain back to time domain." << std::endl;
}
```
convertFrequencyToTimeDomain() utilizes IFFT1D to revert m_fftResult, 
representing frequency domain data, back into time domain PCM samples (m_pcmSamples). 
This process reconstructs the original audio waveform from its frequency domain form,
enabling manipulations such as tone removal in the frequency domain before conversion to the time domain for playback or further processing.
#

6. Visualization: Renders PCM waveform and frequency domain representation using SDL.

```cpp
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
```
clears the renderer, sets the draw color to green, and draws lines connecting PCM samples scaled to fit within a window. 
This visualization helps users see the audio waveform graphically and is one of the requirements.

```cpp
void PCMVisualizer::renderFrequencyDomain() {
SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
SDL_RenderClear(m_renderer);

    // Random color for frequency domain lines
    Uint8 r = rand() % 256;
    Uint8 g = rand() % 256;
    Uint8 b = rand() % 256;
    SDL_SetRenderDrawColor(m_renderer, r, g, b, 255);

    int mid = WINDOW_HEIGHT / 2;
    size_t numBins = m_fftResult.size() / 2;
    float x_scale = static_cast<float>(WINDOW_WIDTH) / static_cast<float>(numBins);

    float maxMagnitude = 0.0f;
    size_t maxIndex = 0;
    for (size_t i = 0; i < numBins; ++i) {
        float magnitude = std::abs(m_fftResult[i]);
        if (magnitude > maxMagnitude) {
            maxMagnitude = magnitude;
            maxIndex = i;
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

    double binWidth = m_sampleRate / static_cast<double>(m_fftResult.size());
    double peakFrequency = maxIndex * binWidth;

    std::cout << "Largest spike frequency: " << peakFrequency << " Hz" << std::endl;
}
```
This visualizes the frequency components of audio data using SDL.
It clears the renderer, draws colored lines representing FFT results scaled to fit the window,
and identifies the highest peak frequency in the spectrum for analysis.

7. Output: Writes modified PCM data to a new file called "AUDIO-AFTER.pcm".

```cpp
void PCMVisualizer::writePCMToFile() {
    std::ofstream outFile("AUDIO-AFTER.pcm", std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening output file: AUDIO-AFTER.pcm" << std::endl;
        return;
    }

    // Write PCM samples to the file
    for (size_t i = 0; i < m_pcmSamples.size(); ++i) {
        outFile.write(reinterpret_cast<const char*>(&m_pcmSamples[i]), sizeof(short));
    }

    outFile.close();
}
```

### Aside

If the user was processing a file with multiple frequencies of interferance,
they could run the code multiple times using the file being created as the next runs input.

The Readme has been formatted with Markdown because I have been using github to save my progress and be able to access from home.

Also the code snippets are included to make understanding it easier and I don't expect them to contribute to word count.