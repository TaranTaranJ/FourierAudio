// Thanks to Imogen Godfrey for help with parsing my array to the DFT1D function.
// Chat GPT was used to help fix a few logic errors

//SDL
#include <SDL2/SDL.h>

// Base
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Maths
#include <random>
#include <math.h>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

using namespace std;

/*
====== Order Of Operation ======

1) convert the img.txt file to an array of complex numbers
2) Convert the array into an image
3) Display the image
4) Save the image

5) Transform each row of the complex array using Fourier Transform
6) Show the transformation as an image
7) Fourier Transform back into original image

...
order 66) Pain... order was lost near the end...
================================
*/

// ====== Initialize Structs ======
struct complexArr
{
    float real;
    float imag;
}; 
// struct name had to be changed from "complex", the name I originally gave it because after many hours of troubleshooting
// I found that it was being called in some unrelated VSCode plugin, which I WAS not notified of...

// window size should be the same as the image being processed
constexpr int WinWidth = 512;
constexpr int WinHeight = 512;

complexArr picture[WinWidth][WinHeight];


/*====================================================================================================================*/


// Display the array/draw image (Operation 2 & 3)
void drawScene(SDL_Renderer *renderer, int WinWidth, int WinHeight, int x, int y, int lightColor)
{
    SDL_SetRenderDrawColor(renderer, lightColor, lightColor, lightColor, 255);
    SDL_RenderDrawPoint(renderer, x, y);
}

// Function to save the generated image (Operation 4)
bool saveImageToFile(SDL_Surface* surface, const std::string& filename)
{
    if (SDL_SaveBMP(surface, filename.c_str()) != 0)
    {
        cerr << "SDL save BMP failed: " << SDL_GetError() << endl;
        return false;
    }

    cout << "Image saved to " << filename << endl;
    return true;
}   

// made with help from imogen to be able to parse complex arrays 
// (the amount of methods for storing and parsing a complex array's that i've tried are ludicrous)
// reads the contents of the .txt file line by line, parsing the numeric values and storing them into the 2D vector of "complexArr"
std::vector<std::vector<complexArr>> readNumbers(const std::string& filename) 
{
    std::ifstream file(filename);

    if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return {};
    }

    std::vector<std::vector<complexArr>> numbers2D;
    std::vector<complexArr> currentRow;

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream stream(line);
        stream >> std::ws;

        float number;
        while (stream >> number) {  // Read multiple numbers per line
            currentRow.push_back({number});
            if (currentRow.size() == WinWidth) {  // Check for 512 items in a row
                numbers2D.push_back(currentRow);  // Add the row to the 2D vector
                currentRow.clear();  // Start a new row
            }
        }
    }

    // Add the last row if it's not empty
    if (!currentRow.empty()) {
        numbers2D.push_back(currentRow);
    }
    file.close();
    return numbers2D;
}

// Function to perform Fourier Transform on each row of the complex array (Operation 5)
// This function was made from adapting psuedo code from scratch pixel
// https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/fourier-transform/fourier-transform-intro.html
void DFT1D(const int N, const complexArr *in, complexArr *out)
{ // y = k & x = n (this was simplified for my understanding & reduce upper/ lowercase conflict)
    for (int y = 0; y < N; ++y)
    {
        out[y].real = out[y].imag = 0;// Initialize the real and imaginary parts

        for (int x = 0; x < N; ++x)
        {
            double angle = 2 * M_PI * y * x / N;  // Correct order of x and y
            out[y].real += in[x].real * ( cos(angle)) - in[x].imag * ( sin(angle)); // Real part
            out[y].imag += in[x].real * (-sin(angle)) + in[x].imag * ( cos(angle)); // Imaginary part
        }
    }
}

// Function to perform Inverse Fourier Transform on each row of the complex array (Operation 7)
// This function was made from adapting psuedo code from scratch pixel
// https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/fourier-transform/fourier-transform-intro.html
void iDFT1D(const int N, const complexArr *in, complexArr *out)
{ // y = k & x = n (this was simplified for my understanding/ sanity & reduce upper/ lowercase conflict)
    for (int y = 0; y < N; ++y) 
    {
        out[y].real = 0, out[y].imag = 0;
        for (int x = 0; x < N; ++x) 
        {
            double angle = 2 * M_PI * y * x / N; // Correct order of x and y
            out[y].real += in[x].real * (cos(angle)) - in[x].imag * (sin(angle)); // Real part
            out[y].imag += in[x].real * (sin(angle)) + in[x].imag * (cos(angle)); // Imaginary part
        }
        out[y].real /= N;
        out[y].imag /= N;
    }
}

auto main() -> int
{
    SDL_Event event;
    SDL_Renderer* renderer;
    SDL_Window* window;
    int x;
    int y;

    // Create an SDL_Surface for rendering
    SDL_Surface* surface = SDL_CreateRGBSurface(0, WinWidth, WinHeight, 32, 0, 0, 0, 0);

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(WinWidth, WinHeight, 0, &window, &renderer);
    SDL_SetWindowTitle(window, "RayTraces");


    // ============ Load your image into the Array (Operation 1) ============

    // ====== Transform each row of the complex array using Fourier Transform (Operation 5) ======
    int N = WinWidth;
    std::string filename = "img.txt";
    // call the readNumbers function
    std::vector<std::vector<complexArr>> numbers = readNumbers(filename);
    
    complexArr in[WinHeight][WinWidth];
    complexArr out[WinHeight][WinWidth];

    for (int i = 0; i < WinHeight; i++) 
    {
        for (int j = 0; j < WinWidth; j++) 
        {
            in[i][j].real = numbers[i][j].real;
        }
    }

    for (int y = 0; y < WinHeight; ++y)
    {
        // Perform Fourier Transform on each row of the complex array
        DFT1D(N, in[y], out[y]);
    }

    // ======================= Display the array/draw image (Operation 2 & 3) =======================

    // Display original corrupted image
    for (int y = 0; y < WinHeight; y++) {
        for (int x = 0; x < WinWidth; x++) {
            float lightColor = in[y][x].real;
            
            drawScene(renderer, WinWidth, WinHeight, x, y, lightColor);
        }
    }
    // Save image of the base image
    SDL_RenderReadPixels(renderer, NULL, surface->format->format, surface->pixels, surface->pitch);
    SDL_SaveBMP(surface, "base.bmp");
    SDL_FreeSurface(surface);
    
    cout << "Render call 1" << endl;
    // flip buffers
    SDL_RenderPresent(renderer);
    SDL_Delay(2000);

    // Render using the values
    for (int y = 0; y < WinHeight; y++) {
        for (int x = 0; x < WinWidth; x++) {
            float lightColor = out[y][x].real;
            
            drawScene(renderer, WinWidth, WinHeight, x, y, lightColor);
        }
    }
    cout << "Render call 2" << endl;
    // flip buffers
    SDL_RenderPresent(renderer);
    SDL_Delay(2000);

    // InverseFourier Function
    for (int y = 0; y < WinHeight; ++y) 
    {
        iDFT1D(WinHeight, in[y], out[y]);
    }

    // Render using the values
    for (int y = 0; y < WinHeight; y++) {
        for (int x = 0; x < WinWidth; x++) {
            float lightColor = out[y][x].real;
            
            drawScene(renderer, WinWidth, WinHeight, x, y, lightColor);
        }
    }

    cout << "Render call 3" << endl;
    // flip buffers
    SDL_RenderPresent(renderer);
    SDL_Delay(2000);

    // Save a screenshot (Operation 4)
    SDL_RenderReadPixels(renderer, NULL, surface->format->format, surface->pixels, surface->pitch);
    SDL_SaveBMP(surface, "screenshot.bmp");
    SDL_FreeSurface(surface);

    // Main loop
    bool quit = false;
    while (!quit)
    {
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                quit = true;
            }
            else if (event.type == SDL_KEYDOWN)
            {
                // Exits if ESC is pressed
                if (event.key.keysym.sym == SDLK_ESCAPE)
                {
                    quit = true;
                }
                // Saves the image if "s" pressed
                else if (event.key.keysym.sym == SDLK_s)
                {
                    SDL_Surface* surface = SDL_CreateRGBSurface(0, WinWidth, WinHeight, 32, 0, 0, 0, 0);
                    SDL_RenderReadPixels(renderer, NULL, surface->format->format, surface->pixels, surface->pitch);
                    SDL_SaveBMP(surface, "squirrel.bmp");
                    SDL_FreeSurface(surface);
                }
            }
        }
        // flip buffers
        SDL_RenderPresent(renderer);
        // wait 1 ms
        SDL_Delay(1);
    }

    // clean up
    SDL_FreeSurface(surface);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return EXIT_SUCCESS;
}