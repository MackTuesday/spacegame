// Ideas
// Sublight engines, costly to change directions moving between stars.
//  So just go slow. What difference does it make? If we do it this way,
//  there needs to be a reason to go as quickly as possible between stars.
//   Limited character lifespan?
//   Limited environment lifespan?
//   But these ideas might take away the fun of exploring.
// Sublight engines, costly to accelerate, so it's important to make sure
//  there's fuel at the destination. Or store fuel in case there is none.


//#include <windows.h>
#include "math.h"
#include "stdio.h"
#include "SDL.h"
#include "lodepng.h"
#include "SDLCharGraphics.h"
#include "BrentsRand.h"
#include "GlobalConstants.h"
#include "StarSystem.h"
#include <algorithm>
#include <vector>


#define kBitsPerStellarSquare 3
#define kISMapColumns        25
#define kISMapRows           25
#define kIPMapColumns        25
#define kIPMapRows           25
#define kPlanetMapColumns   500
#define kPlanetMapRows      500
#define kInfoColumns         27
#define kInfoRows            25

// These need to be >= any of the other character Columns and Rows #defines.
// Disregard the defines that refer to number of pixel columsn and rows
#define kMaxColumns          50
#define kMaxRows             50


static constexpr uint64_t gOrigin = uint64_t(1) << 63;


struct GameState {
    enum {
        kInterstellar,
        kInterplanetary,
        kPlanetary,
        kNumGameStates,
        kNullGameState
    };
    GameState() : mMapPane(nullptr), mInfoPane(nullptr), mStarSystem(nullptr),
                  mPlanet(nullptr), mX(0), mY(0), mID(kNullGameState)  {}
    SDLCharGraphicsPane* mMapPane;
    SDLCharGraphicsPane* mInfoPane;
    StarSystem* mStarSystem;
    Planet* mPlanet;
    uint64_t mX;
    uint64_t mY;
    uint64_t mStepSize;
    int mLevel;
    int mLastLevel;
    GameState* mLastState;
    unsigned mID;
};


static StarSystem* gCurrentSystem = nullptr;
static Planet* gCurrentPlanet = nullptr;
static StarSystem* gSystemUnderReticle = nullptr;
static GameState* gInterstellarState = new GameState;
static GameState* gInterplanetaryState = new GameState;
static GameState* gPlanetaryState = new GameState;
static GameState* gCurrentState = nullptr;
static unsigned cgBuffer[kMaxColumns*kMaxRows];
static unsigned fgClrBuffer[kMaxColumns*kMaxRows];
static unsigned bgClrBuffer[kMaxColumns*kMaxRows];
static float gTerrainBuffer[kPlanetMapColumns*kPlanetMapRows];
static double gWaterline = 0.0;


//GameState* 
void HandleKeyDown(SDL_Event& sdlEvent, SDL_Window* win, SDLCharGraphics& charGraphics);//,
                         //GameState* state);
void ClearCharGraphicsBuffer(unsigned* buffer, unsigned bufferSize);
void FillColorBuffer(unsigned* buffer, unsigned bufferSize, unsigned color);
void CopyStrToIntBuffer(unsigned* buffer, char* str, unsigned maxlen);
unsigned TempToRGB(double temp, double extraSaturation);
unsigned LumToSymbol(double lum);
unsigned PlanetMassToSymbol(double massInKg);
unsigned LifeRGB(double life);
unsigned PlanetSymbolRGB(Planet* planet, StarSystem* starSystem);
//vector<unsigned> PlanetPalette(uint64_t seed, unsigned size);
void FormatInfo(unsigned* buffer, uint64_t posX, uint64_t posY, bool star, 
                uint64_t starTableMod, uint64_t starTableLookup);
void PaintTerrain(SDL_Window* w);
int UpdateGraphics(SDLCharGraphics& cg, SDL_Window* win);//, GameState* state);
inline double PerlinNoise(double x, double y);
inline double PerlinNoise(double x, double y, double z);
inline double Hermite(double x);
inline uint64_t ToFixedPoint12p52(double x);
inline uint64_t NoiseFunction1(uint64_t x);
inline uint64_t NoiseFunction2(uint64_t x);


int main(int argc, char** argv)
{
    //AllocConsole();
    //HANDLE theConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    //freopen("CON", "w", stdout);

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* win = SDL_CreateWindow("roguelike", 40, 40, 1000, 500, 0);

    std::vector<unsigned char> fontSheetVec, fontSheetVec2, fontSheetVec3;
    unsigned fontSheetWidth, fontSheetHeight;
    lodepng::decode(fontSheetVec, fontSheetWidth, fontSheetHeight, "fontsheet.png");
    unsigned* fontSheetPixels = (unsigned*)fontSheetVec.data();
    for (unsigned i = 0; i < fontSheetVec.size()/4; i++) {
        unsigned val = fontSheetPixels[i];
        unsigned swapped = ((val & 0xff000000) >> 24) | ((val & 0x00ff0000) >> 8) |
                           ((val & 0x0000ff00) << 8)  | ((val & 0x000000ff) << 24);
        fontSheetPixels[i] = swapped;
    }
    lodepng::decode(fontSheetVec2, fontSheetWidth, fontSheetHeight, "fontsheet2.png");
    unsigned* fontSheetPixels2 = (unsigned*)fontSheetVec2.data();
    for (unsigned i = 0; i < fontSheetVec2.size()/4; i++) {
        unsigned val = fontSheetPixels2[i];
        unsigned swapped = ((val & 0xff000000) >> 24) | ((val & 0x00ff0000) >> 8) |
                           ((val & 0x0000ff00) << 8)  | ((val & 0x000000ff) << 24);
        fontSheetPixels2[i] = swapped;
    }
/*  lodepng::decode(fontSheetVec3, fontSheetWidth, fontSheetHeight, "fontsheet3.png");
    unsigned* fontSheetPixels3 = (unsigned*)fontSheetVec3.data();
    for (unsigned i = 0; i < fontSheetVec3.size()/4; i++) {
        unsigned val = fontSheetPixels3[i];
        unsigned swapped = ((val & 0xff000000) >> 24) | ((val & 0x00ff0000) >> 8) |
                           ((val & 0x0000ff00) << 8)  | ((val & 0x000000ff) << 24);
        fontSheetPixels3[i] = swapped;
    } */

    SDLCharGraphics charGraphics;

    gInterstellarState->mMapPane = 
        new SDLCharGraphicsPane(kISMapColumns, kISMapRows, 16, 16, 2, 2, 
                                fontSheetPixels, 16, 22, 0, 0, true);
    gInterstellarState->mInfoPane =
        new SDLCharGraphicsPane(kInfoColumns, kInfoRows, 14, 16, 2, 2, 
                                fontSheetPixels2, 16, 22, 500, 0, false);
    gInterstellarState->mStarSystem = nullptr;
    gInterstellarState->mPlanet = nullptr;
    gInterstellarState->mX = 0;
    gInterstellarState->mY = 0;
    gInterstellarState->mStepSize = uint64_t(1) << (32-kBitsPerStellarSquare);
    gInterstellarState->mLevel = 0;
    gInterstellarState->mLastLevel = 0;
    gInterstellarState->mLastState = nullptr;
    gInterstellarState->mID = GameState::kInterstellar;

    gInterplanetaryState->mMapPane = 
        new SDLCharGraphicsPane(kIPMapColumns, kIPMapRows, 16, 16, 2, 2, 
                                fontSheetPixels, 16, 22, 0, 0, true);
    gInterplanetaryState->mInfoPane =
        new SDLCharGraphicsPane(kInfoColumns, kInfoRows, 14, 16, 2, 2, 
                                fontSheetPixels2, 16, 22, 500, 0, false);
    gInterplanetaryState->mStarSystem = nullptr;
    gInterplanetaryState->mPlanet = nullptr;
    gInterplanetaryState->mX = 0;
    gInterplanetaryState->mY = 0;
    gInterplanetaryState->mStepSize = 1;
    gInterplanetaryState->mLevel = 0;
    gInterplanetaryState->mLastLevel = 0;
    gInterplanetaryState->mLastState = nullptr;
    gInterplanetaryState->mID = GameState::kInterplanetary;

    gPlanetaryState->mMapPane = nullptr;
    gPlanetaryState->mInfoPane =
        new SDLCharGraphicsPane(kInfoColumns, kInfoRows, 14, 16, 2, 2, 
                                fontSheetPixels2, 16, 22, 500, 0, false);
    gPlanetaryState->mStarSystem = nullptr;
    gPlanetaryState->mPlanet = nullptr;
    gPlanetaryState->mX = 0;
    gPlanetaryState->mY = 0;
    gPlanetaryState->mStepSize = uint64_t(1) << 21; // Max 2^26 meters, divided into 2^5 steps
    gPlanetaryState->mLevel = 0;
    gPlanetaryState->mLastLevel = 0;
    gPlanetaryState->mLastState = nullptr;
    gPlanetaryState->mID = GameState::kPlanetary;
    
    gCurrentState = gInterstellarState;

    gInterstellarState->mX = uint64_t(1312) << 25;
    gInterstellarState->mY = (uint64_t(549755731952) << 25) + (16 << 25);
    UpdateGraphics(charGraphics, win);//, gCurrentState);

    SDL_Event sdlEvent;
    do {
        gCurrentState->mLastState = gCurrentState;
        if (0 != SDL_PollEvent(&sdlEvent)) {
            if (SDL_KEYDOWN == sdlEvent.type) {
                /*gCurrentState =*/ HandleKeyDown(sdlEvent, win, charGraphics);//, gCurrentState);
            }
            else if (SDL_WINDOWEVENT == sdlEvent.type &&
                     (SDL_WINDOWEVENT_EXPOSED == sdlEvent.window.event ||
                      SDL_WINDOWEVENT_MAXIMIZED == sdlEvent.window.event ||
                      SDL_WINDOWEVENT_RESTORED == sdlEvent.window.event ||
                      SDL_WINDOWEVENT_FOCUS_GAINED == sdlEvent.window.event)) {
                UpdateGraphics(charGraphics, win);//, gCurrentState);
            }
        }
    }
    while (gCurrentState != nullptr);
    
    delete gInterstellarState->mMapPane;
    delete gInterstellarState->mInfoPane;
    delete gInterplanetaryState->mMapPane;
    delete gInterplanetaryState->mInfoPane;
    delete gPlanetaryState->mInfoPane;

    SDL_Quit();
    //FreeConsole();

    return 0;
}


//GameState* 
void HandleKeyDown(SDL_Event& sdlEvent, SDL_Window* win, SDLCharGraphics& charGraphics)//,
                         //GameState* state)
{
    uint64_t d = (gCurrentState->mLevel > 0) ?
                 gCurrentState->mStepSize <<  gCurrentState->mLevel :
                 gCurrentState->mStepSize >> -gCurrentState->mLevel;
    gCurrentState->mLastLevel = gCurrentState->mLevel;
    GameState* nextState = gCurrentState;
    bool updateGraphics = true;
    switch (sdlEvent.key.keysym.sym) {
        case SDLK_KP_1:
        gCurrentState->mX -= d;
        gCurrentState->mY -= d;
        break;

        case SDLK_KP_2:
        case SDLK_DOWN:
        gCurrentState->mY -= d;
        break;

        case SDLK_KP_3:
        gCurrentState->mX += d;
        gCurrentState->mY -= d;
        break;

        case SDLK_KP_4:
        case SDLK_LEFT:
        gCurrentState->mX -= d;
        break;

        case SDLK_KP_6:
        case SDLK_RIGHT:
        gCurrentState->mX += d;
        break;

        case SDLK_KP_7:
        gCurrentState->mX -= d;
        gCurrentState->mY += d;
        break;

        case SDLK_KP_8:
        case SDLK_UP:
        gCurrentState->mY += d;
        break;

        case SDLK_KP_9:
        gCurrentState->mX += d;
        gCurrentState->mY += d;
        break;
            
        case SDLK_PERIOD:
        if (sdlEvent.key.keysym.mod & (KMOD_LSHIFT | KMOD_RSHIFT)) {
            if (GameState::kInterstellar == gCurrentState->mID && 
                gSystemUnderReticle != nullptr) {
                nextState = gInterplanetaryState;
                nextState->mLastState = gCurrentState;
            }
            else if (GameState::kInterplanetary == gCurrentState->mID) {
                nextState = gPlanetaryState;
                nextState->mLastState = gCurrentState;
                nextState->mX = 0;
                nextState->mY = 0;
            }
            else {
                gCurrentState->mLevel--;
            }
        }
        break;
        
        case SDLK_COMMA:
        if (sdlEvent.key.keysym.mod & (KMOD_LSHIFT | KMOD_RSHIFT)) {
            if (GameState::kPlanetary == gCurrentState->mID) {
                if (0 == gCurrentState->mLevel) {
                    nextState = gInterplanetaryState;
                    nextState->mLastState = gCurrentState;
                } else {
                    gCurrentState->mLevel++;
                }
            }
            else if (GameState::kInterplanetary == gCurrentState->mID) {
                nextState = gInterstellarState;
                nextState->mLastState = gCurrentState;
            }
        }
        break;

        default:
        updateGraphics = false;
        break;
    }
    
    if (updateGraphics) {
        gCurrentState = nextState;
        UpdateGraphics(charGraphics, win);//, nextState);
    }
    gCurrentState = (SDLK_END == sdlEvent.key.keysym.sym) ? nullptr : gCurrentState;
    //return nextState;
}
    

void ClearCharGraphicsBuffer(unsigned* buffer, unsigned bufferSize)
{
    for (unsigned i = 0; i < bufferSize; i++) {
        buffer[i] = 32;
    }
}


void FillColorBuffer(unsigned* buffer, unsigned bufferSize, unsigned color)
{
    for (unsigned i = 0; i < bufferSize; i++) {
        buffer[i] = color;
    }
}


void CopyStrToIntBuffer(unsigned* buffer, char* str, unsigned maxlen)
{
    unsigned index = 0;
    while (str[index] != '\0' && index < maxlen) {
        buffer[index] = str[index];
        index++;
    }
}


unsigned TempToRGB(double temp, double extraSaturation)
{
    //RGB components
    double red = 0.0;
    double grn = 0.0;
    double blu = 0.0;
    
    double tempOver100 = temp * 0.01;

    if (temp <= 6600.0) {
        red = 255.0;
    }
    else {
        red = 329.698727446 * pow(tempOver100 - 60.0, -0.1332047592);
    }

    if (red < 0.0) {
        red = 0.0;
    }
    if (red > 255.0) {
        red = 255.0;
    }

    if (temp <= 6600.0) {
        grn = 99.4708025861 * log(tempOver100) - 161.1195681661;
    }
    else {
	    grn = 288.1221695283 * pow(tempOver100 - 60.0, -0.0755148492);
    }

    if (grn < 0.0) {
        grn = 0.0;
    }
    if (grn > 255.0) {
        grn = 255.0;
    }

    if (tempOver100 >= 66.0) {
        blu = 255.0;
    }
    else {
	    if (tempOver100 <= 19.0) {
            blu = 0.0;
        }
        else {
		    blu = 138.5177312231 * log(tempOver100 - 10.0) - 305.0447927307;
	    }
    }

    if (blu < 0.0) {
        blu = 0.0;
    }
    if (blu > 255.0) {
        blu = 255.0;
    }

    //Kelvin Temperature to RGB | Credit to Tanner Helland for the base algorithm
    //Port to C++ for Unreal Engine by Jorge Valle Hurtado - byValle
    //Adapted by me 31 Jan 22
    
    if (extraSaturation > 0.0) {
        double minChannel = red;
        if (grn < minChannel) {
            minChannel = grn;
        }
        if (blu < minChannel) {
            minChannel = blu;
        }
        if (red == minChannel) {
            red /= extraSaturation + 1.0;
        }
        if (grn == minChannel) {
            grn /= extraSaturation + 1.0;
        }
        if (blu == minChannel) {
            blu /= extraSaturation + 1.0;
        }
    }
    
    unsigned r = unsigned(floor(red + 0.5));
    unsigned g = unsigned(floor(grn + 0.5));
    unsigned b = unsigned(floor(blu + 0.5));
    unsigned rgb = (r << 16) | (g << 8) | b;
    
    return rgb;
}


unsigned LumToSymbol(double lum)
{
    unsigned result = 42;
    if (lum < 0.08) {
        result = 13;
    }
    else if (lum < 2.8) {
        result = 12;
    }
    else if (lum < 25) {
        result = 43;
    }
    return result;
}


unsigned PlanetMassToSymbol(double massInKg) {
    // We expect log10 in [23.27, 28.39]
    unsigned symbols[8] {13, 22, 12, 23, 24, 25, 26, 27};
    if (massInKg <= 0.0) {
        return symbols[0];
    }
    constexpr double range = 28.39 - 23.27;
    int index = int(floor((log(massInKg) / log(10.0) - 23.27) / range * 8.0));
    index = (index >= 0) ? index : 0;
    index = (index <= 7) ? index : 7;
    return symbols[index];
}


unsigned LifeRGB(double life)
{
    double hue = fmod(life * 1.0e6, 3.0);
    unsigned r = 0;
    unsigned g = 0;
    unsigned b = 0;
    double f = fmod(hue, 1.0);
    unsigned x = unsigned(floor(256.0 * f));
    unsigned y = 255 - x;
    if (hue < 0.5) {
        r = 255;
        g = x;
        b = 0;
    }
    else if (hue < 1.0) {
        r = y;
        g = 255;
        b = 0;
    }
    else if (hue < 1.5) {
        r = 0;
        g = 255;
        b = x;
    }
    else if (hue < 2.0) {
        r = 0;
        g = y;
        b = 255;
    }
    else if (hue < 2.5) {
        r = x;
        g = 0;
        b = 255;
    }
    else {
        r = 255;
        g = 0;
        b = y;
    }
    return (r << 16) | (g << 8) | b;
}


unsigned PlanetSymbolRGB(Planet* planet, StarSystem* starSystem)
{
    if (planet->mLife >= 2.0) {
        return LifeRGB(planet->mLife);
    }
    else if (planet->mPressure > 1000.0) {
        if (planet->mTemperature > 100) {
            return 0x00ffc38d;
        }
        else {
            return 0x004dbeff;
        }
    }
    else if (planet->mCloudCover > 0.5) {
        if (planet->mTemperature > 500) {
            return 0x00fff5bc;
        }
        else if (planet->mHydCover > 0.5) {
            return 0x00a2a2ff;
        }
        else {
            return 0x00ffc487;
        }
    }
    else if (planet->mTemperature > 1000) {
        return 0x007f0000;
    }
    else if (planet->mFrostCover > 0.5) {
        return 0x00bfbfbf;
    }
    else if (planet->mHydCover > 0.5) {
        return 0x003636ff;
    }
    else {
        return 0x00845c18;
    }
}

/*
vector<unsigned> PlanetPalette(uint64_t seed, unsigned size)
{
    vector<unsigned> result;
    unsigned i = 0;
    for ( ; i < size; i+=16) {
        
    }
    return result;
}
*/

void FormatInfo(unsigned* buffer, uint64_t posX, uint64_t posY, StarSystem* starSystem)
{
    char line[32];
    sprintf(line, "%012lu x %012lu", posX >> 25, posY >> 25);
    CopyStrToIntBuffer(buffer, line, kInfoColumns);
    
    if (starSystem != nullptr) {
        Star* star = starSystem->mStar;
        sprintf(line, "%s%s %5.0lf %7.1le", 
                star->mSpectralClass.c_str(), star->mLuminosityClass.c_str(),
                star->mTemperature, star->mLuminosity);
        CopyStrToIntBuffer(buffer+kInfoColumns, line, kInfoColumns);
        unsigned row = 2;
        for (auto p: starSystem->mPlanets) {
            sprintf(line, "%5.2lf %7.2lf %1.1lf", 
                    p->mSemimajorAxis, p->mMass/kMassOfEarth, p->mLife);
            CopyStrToIntBuffer(buffer+kInfoColumns*row, line, kInfoColumns);
            row++;
            if (row >= kInfoRows) {
                break;
            }
        }
    }
}


void PaintTerrain(SDL_Window* w)
{
    int terrainCols = 500;
    int terrainRows = 250;//(gCurrentState->mLevel == 0) ? 250 : 500;
    if (gPlanetaryState->mLastState != gPlanetaryState ||
        gPlanetaryState->mLastLevel != gPlanetaryState->mLevel) {
        uint64_t paramSeed = 
            NoiseFunction1(uint64_t(fmod(gCurrentPlanet->mSemimajorAxis, 1.0) * 4294967296.0));
//printf("%lu ", paramSeed);
        double a = fmod(double(paramSeed >> 13) / 16777216.0, 1.0);
//printf"%lf ", double(paramSeed >> 13) / 16777216.0);
        paramSeed = NoiseFunction1(paramSeed);
        double b = fmod(double(paramSeed >> 13) / 16777216.0, 1.0);
        paramSeed = NoiseFunction1(paramSeed);
        // If c or d < 0.5 or so, the terrain takes on a somewhat artificial quality.
        double c = fmod(double(paramSeed >> 13) / 16777216.0, 1.0) * 4.0 + 0.5;
        paramSeed = NoiseFunction1(paramSeed);
        double d = fmod(double(paramSeed >> 13) / 16777216.0, 1.0) * 4.0 + 0.5;
        double curveAtZero = -b * pow(a, c) + (1.0-b) * pow(a, d);
        double curveMinimum = 0.0;
        // Ignore the optimum if it occurs at x < 0
        if ((-pow((b*c / ((1.0-b)*d)), 1.0/(d-c)) + a) / (a+1.0) >= 0.0) {
            double curveAtOptimum = pow(b*c / ((1.0-b)*d), c/(d-c)) * (-b + b*c/d);
            curveMinimum = min(curveMinimum, curveAtOptimum);
        }
        curveMinimum = min(curveMinimum, curveAtZero);
printf("%.4lf %.4lf %.4lf %.4lf %.4lf\n", a, b, c, d, curveMinimum);
        double curveNormalizer = 1.0 / (1.0 - curveMinimum);
        paramSeed = NoiseFunction1(paramSeed);
        // 67108864 = 2^26 and 26 = 52/2 and 52 is the number of bits in the
        // mantissa of double precision floating point. Also 52+12 = 64.
        double planetx = double(paramSeed >> 12) / 67108864.0 + 67108865.0;
        paramSeed = NoiseFunction1(paramSeed);
        double planety = double(paramSeed >> 12) / 67108864.0 + 67108865.0;
        paramSeed = NoiseFunction1(paramSeed);
        double planetz = double(paramSeed >> 12) / 67108864.0 + 67108865.0;
        paramSeed = NoiseFunction1(paramSeed);
//printf("%.4lf %.4lf %.4lf\n", planetx, planety, planetz);
        vector<float> analysisBuffer;
        unsigned numLayers = 10;
        double coordScales[20];
        coordScales[0] = 0.3;
        for (unsigned n = 1; n < numLayers; ++n) {
            coordScales[n] = coordScales[n-1] *
                             (fmod(coordScales[n-1]*1.6473489695, 2.2) + 1.1);
        }
        double warpScale = 0.75;
        double warpStrength = fmod(double(paramSeed >> 12) / 65536.0, 1.0) * 2.5;
        for (int j = 0; j < terrainRows; ++j) {
            double lat = 
                kPI * (0.5 - double(j) / 250 + 
                       (gCurrentState->mY & 31) / 32.0) /
                (1 << -gCurrentState->mLevel);
            for (int i = 0; i < terrainCols; ++i) {
                unsigned index = j*terrainCols + i;
                double lon = 
                    kPI * (2.0*i/terrainCols - 1.0 + 
                           (gCurrentState->mX & 31) / 16.0) /
                    (1 << -gCurrentState->mLevel);
                double heightScale = 1.0;
                double coslon = cos(lon);
                double sinlon = sin(lon);
                double t = 0.0;
                double x = coslon + planetx - 134217729.0;
                double y = sinlon + planety - 134217729.0;
                double z = lat + planetz - 134217729.0;
                double beta = 0.25 + 0.3 * PerlinNoise(x, y, z);
                double heightNormalizer = (1.0 - beta) / (1.0 - pow(beta, numLayers));
                double warp = (PerlinNoise(warpScale*coslon + planetx + 1000.0,
                                           warpScale*sinlon + planety + 1000.0, 
                                           warpScale*lat + planetz + 1000.0) - 0.5) * 2.0;
                warp *= warp * warpStrength;
                for (unsigned n = 0; n < numLayers; ++n) {
                    double zscale = coordScales[n];
                    x = zscale * coslon;
                    y = zscale * sinlon;
                    z = coordScales[n] * lat;
                    double w = 1.0 + warp / coordScales[n];
                    double h = PerlinNoise(x*w+planetx, y*w+planety, z*w+planetz);
                    double q = h * a + h - a; // [0,1) ==> [-a,1)
                    // [-a,1) ==> [-ba^c + (1-b)a^d, 1)
                    double absq = (q < 0.0) ? -q : q;
                    double logq = log(absq + 1.0e-300);
                    double qToTheD = exp(logq * d);
                    double r = b * (exp(logq*c) * kFsgn(q) - qToTheD) + qToTheD;
                    // [-ba^c + (1-b)a^d, 1) ==> [0,1)
                    double s = (r - curveMinimum) * curveNormalizer;
                    s = (0 == (paramSeed & (0x1 << n))) ? s : 1.0-s;
                    t += s * heightScale;
                    heightScale *= beta;
                }
//printf("%02d %02d %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", i, j, lat, lon, x, y, z, h);
                gTerrainBuffer[index] = float(t * heightNormalizer);
                analysisBuffer.push_back(float(t * heightNormalizer));
            }
        }
        sort(analysisBuffer.begin(), analysisBuffer.end());
        gWaterline = analysisBuffer[unsigned(floor((analysisBuffer.size()-1) *
                                                   gCurrentPlanet->mHydCover))];
    }

    SDL_Surface* surf = SDL_GetWindowSurface(w);
    unsigned* pixelsPtr = ((unsigned*)surf->pixels);

    for (int j = 0; j < 500; ++j) {
        for (int i = 0; i < 500; ++i) {
            pixelsPtr[j * surf->w + i] = 0x00000000;
        }
    }
    
    int signedX = int(gCurrentState->mX & 0xffffffff);
    int signedY = int(gCurrentState->mY & 0xffffffff);
    int sjbeg = (500-terrainRows) / 2;
    int sjend = sjbeg + terrainRows;
    int tjbeg = 0;
    int tjend = terrainRows;
    int tibeg = 0;
    int tiend = terrainCols;
    int sibeg = 0;
    int siend = terrainCols;
    // 500/2^26 = 0.000007450580596923828125
    double moveScale = 7.450580596923828125e-6 * (1 << -gCurrentState->mLevel);
    int cropX = int(floor(signedX * moveScale + 0.5));
    int cropY = int(floor(signedY * moveScale + 0.5));

    if (signedY < 0) {
        tjbeg -= cropY;
        sjend += cropY;
    }
    else {
        sjbeg += cropY;
        tjend -= cropY;
    }
    if (signedX < 0) {
        sibeg -= cropX;
        tiend += cropX;
    }
    else {
        tibeg += cropX;
        siend -= cropX;
    }
    
    for (int j = sjbeg; j < sjend; ++j) {
        for (int i = sibeg; i < siend; ++i) {
            int index = (j-sjbeg+tjbeg) * 500 + i-sibeg+tibeg;
            int adjIndex = (j-sjbeg+tjbeg) * 500 + (i-sibeg+tibeg+1) % 500;
            if (gTerrainBuffer[index] > gWaterline) {
                double diff = floor(5000.0 * (gTerrainBuffer[index] - 
                                              gTerrainBuffer[adjIndex])) + 64.0;
                unsigned s = unsigned((diff < 60.0) ? 32.0 : diff);
                s = (s <= 255) ? s : 255;
                pixelsPtr[j * surf->w + i] = 0x00010101 * s;
            }
            else {
                pixelsPtr[j * surf->w + i] = 64;
            }
        }
    }
}

/*
void PaintTerrain(SDL_Window* w)
{
    static uint64_t lastX = 0;
    static uint64_t lastY = 0;
    
    SDL_Surface* surf = SDL_GetWindowSurface(w);

    unsigned* pixelsPtr = ((unsigned*)surf->pixels);
    if (gCurrentState->mLevel == 0) {
        for (int j = 0; j < 125; ++j) {
            for (int i = 0; i < 500; ++i) {
                pixelsPtr[j * surf->w + i] = 0x00000000;
            }
        }
    }
    // f(x)  = b |q|^c sgn(q) + (1-b) |q|^d
    //       = b |ax-a+x|^c sgn(q) + (1-b) |ax-a+x|^d
    // df/dx = dq/dx (bc|q|^(c-1) d|q|/dq sgn(q) + b|q|^c delta(q) + 
    //                (1-b) d|q|^(d-1) d|q|/dq)
    //       = (a+1) (bc|q|^(c-1) sgn(q) sgn(q) + b|q|^c delta(q) + 
    //                (1-b) d|q|^(d-1) sgn(q))
    //       = (a+1) (bc|q|^(c-1) + (1-b) d|q|^(d-1) sgn(q))             q != 0
    // Find the minimum:
    // 0     = (a+1) (bc|q|^(c-1) + (1-b) d|q|^(d-1) sgn(q))             q != 0
    // 0     = bc + (1-b) d|q|^(d-c) sgn(q)
    // -bc / (1-b)d = |q|^(d-c) sgn(q)
    // The left is negative so sgn(q) must be negative
    // bc / (1-b)d = |q|^(d-c)
    // (bc / (1-b)d)^(1/(d-c)) = |q|
    // ±|-bc / (1-b)d|^(1/(d-c)) = qmin
    //                           = ax-a+x
    // (±|-bc / (1-b)d|^(1/(d-c)) + a) / (a+1) = xmin
    // We want the root with the negative sign because sgn(q) is negative
    // f(qmin) = b  |-|-bc / (1-b)d|^(1/(d-c))|^c sgn(-|-bc / (1-b)d|^(1/(d-c))) + 
    //        (1-b) |-|-bc / (1-b)d|^(1/(d-c))|^d
    //         = -b ||-bc / (1-b)d|^(1/(d-c))|^c + (1-b) ||-bc / (1-b)d|^(1/(d-c))|^d
    //         = -b ||-bc / (1-b)d|^(1/(d-c))|^c + (1-b) ||-bc / (1-b)d|^(1/(d-c))|^d
    //         = -b |-bc / (1-b)d|^(c/(d-c)) + (1-b) |-bc / (1-b)d|^(d/(d-c))
    //         = |-bc / (1-b)d|^(c/(d-c)) 
    //           (-b + (1-b) |-bc / (1-b)d|^(d/(d-c)) |-bc / (1-b)d|^(-c/(d-c)))
    //         = |-bc / (1-b)d|^(c/(d-c)) (-b + (1-b) |-bc / (1-b)d|^((d-c)/(d-c)))
    //         = |-bc / (1-b)d|^(c/(d-c)) (-b + (1-b) |-bc / (1-b)d|)
    // 1-b is non-negative so we can do this
    //         = |-bc / (1-b)d|^(c/(d-c)) (-b + |-bc/d|)
    // b, c, and d are all non-negative, so...
    //         = (bc / (1-b)d)^(c/(d-c)) (-b + bc/d)
    // So we'll compare the values at xmin, q=0 <==> x=(q+a)/(1+a), x=0, x=1.
    // The value at q=0 is 0.
    // f(0) = -ba^c + (1-b) a^d
    // f(1) = 1
    // 1 > 0 so we can disregard value at x=1
    uint64_t paramSeed = 
        NoiseFunction1(uint64_t(fmod(gCurrentPlanet->mSemimajorAxis, 1.0) * 4294967296.0));
//printf("%lu ", paramSeed);
    double a = fmod(double(paramSeed >> 13) / 16777216.0, 1.0);
//printf("%lf ", double(paramSeed >> 13) / 16777216.0);
    paramSeed = NoiseFunction1(paramSeed);
    double b = fmod(double(paramSeed >> 13) / 16777216.0, 1.0);
    paramSeed = NoiseFunction1(paramSeed);
    // If c or d < 0.5 or so, the terrain takes on a somewhat artificial quality.
    double c = fmod(double(paramSeed >> 13) / 16777216.0, 1.0) * 4.0 + 0.5;
    paramSeed = NoiseFunction1(paramSeed);
    double d = fmod(double(paramSeed >> 13) / 16777216.0, 1.0) * 4.0 + 0.5;
    double curveAtZero = -b * pow(a, c) + (1.0-b) * pow(a, d);
    double curveMinimum = 0.0;
    // Ignore the optimum if it occurs at x < 0
    if ((-pow((b*c / ((1.0-b)*d)), 1.0/(d-c)) + a) / (a+1.0) >= 0.0) {
        double curveAtOptimum = pow(b*c / ((1.0-b)*d), c/(d-c)) * (-b + b*c/d);
        curveMinimum = min(curveMinimum, curveAtOptimum);
    }
    curveMinimum = min(curveMinimum, curveAtZero);
//printf("%.4lf %.4lf %.4lf %.4lf %.4lf\n", a, b, c, d, curveMinimum);
    double curveNormalizer = 1.0 / (1.0 - curveMinimum);
    vector<float> analysisBuffer;
    paramSeed = NoiseFunction1(paramSeed);
    // 67108864 = 2^26 and 26 = 52/2 and 52 is the number of bits in the
    // mantissa of double precision floating point. Also 52+12 = 64.
    double planetx = double(paramSeed >> 12) / 67108864.0 + 67108865.0;
    paramSeed = NoiseFunction1(paramSeed);
    double planety = double(paramSeed >> 12) / 67108864.0 + 67108865.0;
    paramSeed = NoiseFunction1(paramSeed);
    double planetz = double(paramSeed >> 12) / 67108864.0 + 67108865.0;
    paramSeed = NoiseFunction1(paramSeed);
//printf("%.4lf %.4lf %.4lf\n", planetx, planety, planetz);
    unsigned numLayers = 10;
    double coordScales[20];
    coordScales[0] = 1.0;
    for (unsigned n = 1; n < numLayers; ++n) {
        coordScales[n] = coordScales[n-1] *
                         (fmod(coordScales[n-1]*1.6473489695, 2.2) + 1.1);
    }
    double warpScale = 0.75;
    double warpStrength = fmod(double(paramSeed >> 12) / 65536.0, 1.0) * 2.5;
    int jbeg = (0 == gCurrentState->mLevel) ? 125 :   0;
    int jend = (0 == gCurrentState->mLevel) ? 376 : 500;
printf("%lu %lu\n", gCurrentState->mX, gCurrentState->mY);
    for (int j = jbeg; j < jend; ++j) {
        double lat = kPI * (0.5 / 126.0 * (250-j) + (gCurrentState->mY & 31) / 32.0) / 
                     (1 << -gCurrentState->mLevel);
        double coslat = cos(lat);
        double sinlat = sin(lat);
        for (int i = 0; i < 500; ++i) {
            unsigned index = (j-jbeg)*500 + i;
            double lon = kPI * (1.0 / 250.0 * (i-250) + (gCurrentState->mX & 31) / 16.0) / 
                         (1 << -gCurrentState->mLevel);
            double coordScale = 1.0;
            double heightScale = 1.0;
            double coslon = cos(lon);
            double sinlon = sin(lon);
            double t = 0.0;
            double x = coslat*coslon + planetx - 134217729.0;
            double y = coslat*sinlon + planety - 134217729.0;
            double z = sinlat + planetz - 134217729.0;
            double beta = 0.25 + 0.3 * PerlinNoise(x, y, z);
            double heightNormalizer = (1.0 - beta) / (1.0 - pow(beta, numLayers));
            double warp = (PerlinNoise(warpScale*coslat*coslon + planetx + 1000.0,
                                       warpScale*coslat*sinlon + planety + 1000.0, 
                                       warpScale*sinlat + planetz + 1000.0) - 0.5) * 2.0;
            warp *= warp * warpStrength;
            for (unsigned n = 0; n < numLayers; ++n) {
                double zscale = coordScales[n] * coslat;
                x = zscale * coslon;
                y = zscale * sinlon;
                z = coordScales[n] * sinlat;
                double w = 1.0 + warp / coordScales[n];
                double h = PerlinNoise(x*w+planetx, y*w+planety, z*w+planetz);
                double q = h * a + h - a; // [0,1) ==> [-a,1)
                // [-a,1) ==> [-ba^c + (1-b)a^d, 1)
                double absq = (q < 0.0) ? -q : q;
                double logq = log(absq + 1.0e-300);
                double qToTheD = exp(logq * d);
                double r = b * (exp(logq*c) * kFsgn(q) - qToTheD) + qToTheD;
                // [-ba^c + (1-b)a^d, 1) ==> [0,1)
                double s = (r - curveMinimum) * curveNormalizer;
                s = (0 == (paramSeed & (0x1 << n))) ? s : 1.0-s;
                t += s * heightScale;
                heightScale *= beta;
            }
//printf("%02d %02d %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", i, j, lat, lon, x, y, z, h);
            gTerrainBuffer[index] = float(t * heightNormalizer);
            analysisBuffer.push_back(float(t * heightNormalizer));
        }
    }
//printf("%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", minh, maxh, minr, maxr, mins, maxs);
    sort(analysisBuffer.begin(), analysisBuffer.end());
    gWaterline = analysisBuffer[unsigned(floor((analysisBuffer.size()-1) *
                                               gCurrentPlanet->mHydCover))];
//printf("%.4lf %.4lf %.4f\n", gCurrentPlanet->mNormedHydDepth, gCurrentPlanet->mHydCover, 
//                             waterline);
    for (int j = jbeg; j < jend; ++j) {
        for (int i = 0; i < 500; ++i) {
            unsigned index = (j-jbeg)*500 + i;
            unsigned adjIndex = (j-jbeg)*500 + (i+1)%500;
            if (gTerrainBuffer[index] > gWaterline) {
                double diff = floor(5000.0 * (gTerrainBuffer[index] - 
                                              gTerrainBuffer[adjIndex])) + 64.0;
                unsigned s = unsigned((diff < 60.0) ? 32.0 : diff);
                s = (s <= 255) ? s : 255;
                pixelsPtr[j * surf->w + i] = 0x00010101 * s;
                //    0x00010101 * unsigned(floor(terrainBuffer[index]*256.0));
            }
            else {
                pixelsPtr[j * surf->w + i] = 64;
            }
        }
    }
    if (gCurrentState->mLevel == 0) {
        for (int j = 376; j < 500; ++j) {
            for (int i = 0; i < 500; ++i) {
                pixelsPtr[j * surf->w + i] = 0x00000000;
            }
        }
    }
}
*/

int UpdateGraphics(SDLCharGraphics& cg, SDL_Window* w)//, GameState* state)
{
printf("UpdateGraphics: state %u\n", gCurrentState->mID);
    SDL_Surface* surf = SDL_GetWindowSurface(w);
    
    // Determine where in space we are
    // Determine contents of that region of space and put them in a buffer for cg
    uint32_t posXInt = (uint32_t)(gCurrentState->mX >> 32);
    uint32_t posXFrac = (uint32_t)(gCurrentState->mX & 0x00000000ffffffff);
    uint32_t posYInt = (uint32_t)(gCurrentState->mY >> 32);
    uint32_t posYFrac = (uint32_t)(gCurrentState->mY & 0x00000000ffffffff);
    if (GameState::kInterstellar == gCurrentState->mID) {
        gSystemUnderReticle = nullptr;
        double step = 1.0 / (double)(1 << kBitsPerStellarSquare);
        for (int j = -12; j < 13; j++) {
            for (int i = -12; i < 13; i++) {
                double x = (double)posXInt + (double)posXFrac/4294967296.0 + i*step;
                double y = (double)posYInt + (double)posYFrac/4294967296.0 + j*step;
                // Adding half a step so the value isn't exactly 0.5 at integer x and y.
                // The amount added mustn't result in any sort of truncation. The sum
                // must be exactly representable in double precision floats even when
                // posXInt and posYInt are large. Otherwise you'll get weird behavior
                // near the wraparound between 2^32 and 0.
                double noiseValue = PerlinNoise(x+step*0.5, y+step*0.5);
                double maxNoiseValue = 1.0 - 1.0/(double)((uint64_t)1 << 53);
                noiseValue = (noiseValue < 0.0) ? 0.0 :
                             ((noiseValue > maxNoiseValue) ? maxNoiseValue : noiseValue);
                double starNoise = noiseValue * ((uint64_t)1 << 43);
                double starFrac = starNoise - floor(starNoise);
                double starFrequency = (63.0 * pow(noiseValue, 8.0) + 1) / 192.0;
                int index = (12-j)*25 + i + 12;
                cgBuffer[index] = 32;
                fgClrBuffer[index] = 0xffffff;
                bgClrBuffer[index] = unsigned(noiseValue * 256.0) | (32 << 16);
                StarSystem* tempSystem = nullptr;
                bool center = ((i == 0) && (j == 0));
                if (starFrac <= starFrequency) {
                    uint64_t fpNoiseValue = ToFixedPoint12p52(noiseValue);
                    tempSystem = new StarSystem(unsigned(fpNoiseValue & 0xffffffff), !center);
                    cgBuffer[index] = LumToSymbol(tempSystem->mStar->mLuminosity);
                    fgClrBuffer[index] = TempToRGB(tempSystem->mStar->mTemperature, 0.5);
                    if (center) {
                        if (gCurrentSystem != nullptr) {
                            delete gCurrentSystem;
                        }
                        gCurrentSystem = tempSystem;
                        gSystemUnderReticle = tempSystem;
                    }
                    else {
                        delete tempSystem;
                    }
                }
            }
        }

        cg.WriteCharacters(gCurrentState->mMapPane, cgBuffer, fgClrBuffer, bgClrBuffer, 0, 0, 
                           kISMapColumns, kISMapRows);
        
        ClearCharGraphicsBuffer(cgBuffer, kMaxColumns*kMaxRows);
        FillColorBuffer(fgClrBuffer, kMaxColumns*kMaxRows, 0xffffff);
        FillColorBuffer(bgClrBuffer, kMaxColumns*kMaxRows, 0);
        FormatInfo(cgBuffer, gCurrentState->mX, gCurrentState->mY, gCurrentSystem);
        cg.WriteCharacters(gCurrentState->mInfoPane, cgBuffer, fgClrBuffer, bgClrBuffer, 0, 0, 
                           kInfoColumns, kInfoRows);
    }
    else if (GameState::kInterplanetary == gCurrentState->mID) {
        ClearCharGraphicsBuffer(cgBuffer, kMaxColumns*kMaxRows);
        FillColorBuffer(fgClrBuffer, kMaxColumns*kMaxRows, 0xffffff);
        FillColorBuffer(bgClrBuffer, kMaxColumns*kMaxRows, 0);
        double minSA = 1e+300;
        double maxSA = -1e+300;
        double minRatio = 1.0e300;
        unsigned numPlanets = gCurrentSystem->mPlanets.size();
        for (unsigned i = 0; i < numPlanets; ++i) {
            Planet* p = gCurrentSystem->mPlanets[i];
            if (p->mSemimajorAxis < minSA) {
                minSA = p->mSemimajorAxis;
            }
            if (p->mSemimajorAxis > maxSA) {
                maxSA = p->mSemimajorAxis;
            }
            if (i > 0) {
                double ratio = p->mSemimajorAxis / 
                               gCurrentSystem->mPlanets[i-1]->mSemimajorAxis;
                if (ratio < minRatio) {
                    minRatio = ratio;
                }
            }
        }
//        double logSARange = 1.0;
        double orbitScale = 1.0;
        if (numPlanets > 1) {
            orbitScale = max(1.01 / log(minRatio), 90.0 / log(maxSA/minSA));
        }
        vector<uint64_t> planetYs;
        uint64_t minY = -1; // It'll underflow to the max unsigned
        for (unsigned i = 0; i < numPlanets; ++i) {
            Planet* p = gCurrentSystem->mPlanets[i];
            uint64_t y = uint64_t(gOrigin - 5) - 
                         uint64_t(floor(orbitScale * (log(p->mSemimajorAxis)-log(minSA))));
            planetYs.push_back(y);
            if (y < minY) {
                minY = y;
            }
//printf("%.4lf %.4lf %.4lf %.4lf\n", p->mSemimajorAxis, p->mPressure, p->mNormedHydDepth, p->mHydCover);
        }
//printf("\n");
        if (gCurrentState->mLastState->mID == GameState::kInterstellar) {
            gCurrentState->mX = gOrigin;
            gCurrentState->mY = minY - 5;
        }
        // Sun
        uint64_t screenX = -gCurrentState->mX + (gOrigin + 12);
        uint64_t screenY =  gCurrentState->mY - (gOrigin - 12);
        if (screenX < 25 && screenY < 25) {
            unsigned index = screenX + 25 * screenY;
            cgBuffer[index] = 42;
            fgClrBuffer[index] = TempToRGB(gCurrentSystem->mStar->mTemperature, 0.5);
            bgClrBuffer[index] = 0x00000000;
        }
        // Planets
        for (unsigned i = 0; i < planetYs.size(); ++i) {
            screenX = -gCurrentState->mX + (gOrigin + 12);
            screenY =  gCurrentState->mY - (planetYs[i] - 12);
            if (screenX < 25 && screenY < 25) {
                unsigned index = screenX + 25 * screenY;
                Planet* tempPlanet = gCurrentSystem->mPlanets[i];
                cgBuffer[index] = PlanetMassToSymbol(tempPlanet->mMass);
                fgClrBuffer[index] = PlanetSymbolRGB(tempPlanet, gCurrentSystem);
                bgClrBuffer[index] = 0x00000000;
                if (screenX == 12 && screenY == 12) {
                    gCurrentPlanet = tempPlanet;
                }
            }
        }
        cg.WriteCharacters(gCurrentState->mMapPane, cgBuffer, fgClrBuffer, bgClrBuffer, 0, 0, 
                           kIPMapColumns, kIPMapRows);
    }
    else {
        PaintTerrain(w);
    }

//    unsigned const* pixelBuffer = cg.GetPixels();
    if (GameState::kPlanetary != gCurrentState->mID) {
        unsigned const* pixelBuffer = cg.GetPixels();
        for (unsigned j = 0; j < 500; j++) {
            for (unsigned i = 0; i < 1000; i++) {
                ((unsigned*)surf->pixels)[j * surf->w + i] = pixelBuffer[j * 1000 + i];
            }
        }
    }

    SDL_UpdateWindowSurface(w);

    return 0;
}


double PerlinNoise(double x, double y)
{
    double xint = floor(x);
    double yint = floor(y);

    uint32_t xuint = (uint32_t)xint;
    uint32_t yuint = (uint32_t)yint;
    
    uint64_t ll = NoiseFunction1(((uint64_t) xuint    << 32) | (uint64_t) yuint);
    uint64_t lr = NoiseFunction1(((uint64_t)(xuint+1) << 32) | (uint64_t) yuint);
    uint64_t ul = NoiseFunction1(((uint64_t) xuint    << 32) | (uint64_t)(yuint+1));
    uint64_t ur = NoiseFunction1(((uint64_t)(xuint+1) << 32) | (uint64_t)(yuint+1));

    // double has 52 bits of mantissa
    double phaseNormalizer = 3.141592653589793 * 2.0 / (double)((uint64_t)1 << 52);
    double llphase = (double)(ll >> (64-52)) * phaseNormalizer;
    double lrphase = (double)(lr >> (64-52)) * phaseNormalizer;
    double ulphase = (double)(ul >> (64-52)) * phaseNormalizer;
    double urphase = (double)(ur >> (64-52)) * phaseNormalizer;

    double xfrac = x - xint;
    double yfrac = y - yint;

    double lldot =      xfrac *cos(llphase) +      yfrac *sin(llphase);
    double lrdot = (xfrac-1.0)*cos(lrphase) +      yfrac *sin(lrphase);
    double uldot =      xfrac *cos(ulphase) + (yfrac-1.0)*sin(ulphase);
    double urdot = (xfrac-1.0)*cos(urphase) + (yfrac-1.0)*sin(urphase);

    double hermitex = Hermite(xfrac);
    double hermitey = Hermite(yfrac);

    double linterp = (1.0-hermitex)*lldot + hermitex*lrdot;
    double uinterp = (1.0-hermitex)*uldot + hermitex*urdot;
    double interp = (1.0-hermitey)*linterp + hermitey*uinterp;

    double result = interp * 0.7071067811865475 + 0.5;
/*
    if (verbose) {
        printf("x = %.16lf, y = %.16lf\n", x, y);
        printf("xint = %.16lf, yint = %.16lf\n", xint, yint);
        printf("xuint = %u, yuint = %u\n", xuint, yuint);
        //printf("ll = %llu, lr = %llu, ul = %llu, ur = %llu\n", ll, lr, ul, ur);
        printf("ll = %lu, lr = %lu, ul = %lu, ur = %lu\n", ll, lr, ul, ur);
        printf("llphase = %.16lf, lrphase = %.16lf, ulphase = %.16lf, urphase = %.16lf\n",
               llphase, lrphase, ulphase, urphase);
        printf("xfrac = %.16lf, yfrac = %.16lf\n", xfrac, yfrac);
        printf("lldot = %.16lf, lrdot = %.16lf, uldot = %.16lf, urdot = %.16lf\n",
               lldot, lrdot, uldot, urdot);
        printf("hermitex = %.16lf, hermitey = %.16lf\n", hermitex, hermitey);
        printf("linterp = %.16lf, uinterp = %.16lf, interp = %.16lf\n",
               linterp, uinterp, interp);
        printf("result = %.16lf\n", result);
        printf("\n");
    }
*/
    return result;
}


inline uint64_t RotRight21(uint64_t n)
{
    return (n >> 21) | (n << 43);
}
inline uint64_t RotRight42(uint64_t n)
{
    return (n >> 42) | (n << 22);
}


inline double GradDotProduct(unsigned n, double x, double y, double z)
{
    switch (n)
    {
        case 0:
        return x + y;
        break;
        
        case 1:
        return y - x;
        break;
        
        case 2:
        return x - y;
        break;
        
        case 3:
        return -x - y;
        break;
      
        case 4:
        return x + z;
        break;
        
        case 5:
        return z - x;
        break;
        
        case 6:
        return x - z;
        break;
        
        case 7:
        return -x - z;
        break;
        
        case 8:
        return y + z;
        break;
        
        case 9:
        return z - y;
        break;
        
        case 10:
        return y - z;
        break;
        
        case 11:
        return -y - z;
        break;
        
        default:
        return 0;
        break;
    }
}


#define kRootHalf 0.7071067811865475
inline double PerlinNoise(double x, double y, double z)
{
    uint64_t xuint = static_cast<uint64_t>(x);
    uint64_t yuint = static_cast<uint64_t>(y);
    uint64_t zuint = static_cast<uint64_t>(z);
    
    uint64_t xuintp1 = xuint + 1;
    uint64_t yuintp1 = yuint + 1;
    uint64_t zuintp1 = zuint + 1;
    
    uint64_t yuintRot21 = RotRight21(yuint);
    uint64_t yuintp1Rot21 = RotRight21(yuintp1);
    
    uint64_t zuintRot42 = RotRight42(zuint);
    uint64_t zuintp1Rot42 = RotRight42(zuintp1);

    unsigned fll = unsigned(NoiseFunction1(xuint   ^ yuintRot21   ^ zuintRot42  ) >> 33);
    unsigned flr = unsigned(NoiseFunction1(xuintp1 ^ yuintRot21   ^ zuintRot42  ) >> 33);
    unsigned ful = unsigned(NoiseFunction1(xuint   ^ yuintp1Rot21 ^ zuintRot42  ) >> 33);
    unsigned fur = unsigned(NoiseFunction1(xuintp1 ^ yuintp1Rot21 ^ zuintRot42  ) >> 33);
    unsigned bll = unsigned(NoiseFunction1(xuint   ^ yuintRot21   ^ zuintp1Rot42) >> 33);
    unsigned blr = unsigned(NoiseFunction1(xuintp1 ^ yuintRot21   ^ zuintp1Rot42) >> 33);
    unsigned bul = unsigned(NoiseFunction1(xuint   ^ yuintp1Rot21 ^ zuintp1Rot42) >> 33);
    unsigned bur = unsigned(NoiseFunction1(xuintp1 ^ yuintp1Rot21 ^ zuintp1Rot42) >> 33);

    // gi stands for gradient index
    unsigned fllgi = (fll + (fll>>1)) >> 28;
    unsigned flrgi = (flr + (flr>>1)) >> 28;
    unsigned fulgi = (ful + (ful>>1)) >> 28;
    unsigned furgi = (fur + (fur>>1)) >> 28;
    unsigned bllgi = (bll + (bll>>1)) >> 28;
    unsigned blrgi = (blr + (blr>>1)) >> 28;
    unsigned bulgi = (bul + (bul>>1)) >> 28;
    unsigned burgi = (bur + (bur>>1)) >> 28;

    double xfrac = x - floor(x);
    double yfrac = y - floor(y);
    double zfrac = z - floor(z);
    
    double xfracm1 = xfrac - 1.0;
    double yfracm1 = yfrac - 1.0;
    double zfracm1 = zfrac - 1.0;

    double flldot = GradDotProduct(fllgi, xfrac,   yfrac,   zfrac);
    double flrdot = GradDotProduct(flrgi, xfracm1, yfrac,   zfrac);
    double fuldot = GradDotProduct(fulgi, xfrac,   yfracm1, zfrac);
    double furdot = GradDotProduct(furgi, xfracm1, yfracm1, zfrac);
    double blldot = GradDotProduct(bllgi, xfrac,   yfrac,   zfracm1);
    double blrdot = GradDotProduct(blrgi, xfracm1, yfrac,   zfracm1);
    double buldot = GradDotProduct(bulgi, xfrac,   yfracm1, zfracm1);
    double burdot = GradDotProduct(burgi, xfracm1, yfracm1, zfracm1);

    double hx = Hermite(xfrac);
    double hy = Hermite(yfrac);
    double hz = Hermite(zfrac);

    double n00x = flldot + hx*(flrdot - flldot);
    double n01x = fuldot + hx*(furdot - fuldot);
    double n10x = blldot + hx*(blrdot - blldot);
    double n11x = buldot + hx*(burdot - buldot);

    double n0yx = n00x + hy*(n01x - n00x);
    double n1yx = n10x + hy*(n11x - n10x);

    double mzyx = n0yx + hz*(n1yx - n0yx);

    double result = mzyx * kRootHalf + 0.5;
    result = (result > 0.0) ? result : 0.0;
    result = (result < 1.0) ? result : 1.0;

    return result;
}


inline uint64_t NoiseFunction1(uint64_t x)
{
    uint64_t result = x * 2862933555777941757 + 3037000493;
    result ^= result << 13;
    result ^= result >> 7;
    result ^= result << 17;
    return result;
}


inline uint64_t NoiseFunction2(uint64_t x)
{
    uint64_t result = x * 6680178296815197433 + 7046029254386353087;
    result ^= result << 13;
    result ^= result >> 7;
    result ^= result << 17;
    return result;
}

// x should be in [0,1]
inline double Hermite(double x)
{
    double xsquared = x * x;
    return x * x * x * (6.0*x * x - 15.0*x + 10.0);
}

// x should be in [0,1]
inline uint64_t ToFixedPoint12p52(double x)
{
    uint64_t crosscast = *((uint64_t*)(&x));
    uint64_t mantissa = 0x0010000000000000 | (crosscast & 0x000fffffffffffff);
    int64_t exponent = (int64_t)((0x7ff0000000000000 & crosscast) >> 52) - 1023;
    exponent = (exponent > 0) ? 0 : exponent;
    exponent = (exponent < -52) ? -52 : exponent;
    uint64_t result = mantissa >> -exponent;
    return result;
}

/*
int UpdateGraphics(SDLCharGraphics& cg, SDL_Window* w, unsigned posX, unsigned posY)
{
    static unsigned charBuffer[625];
    static unsigned fgClrBuffer[625];
    static unsigned bgClrBuffer[625];

    SDL_Surface* s = SDL_GetWindowSurface(w);

    // Determine where in space we are
    // Determine contents of that region of space and put them in a buffer for cg
    for (unsigned j = 0; j < 25; j++) {
        for (unsigned i = 0; i < 25; i++) {
            uint64_t noiseValue = noiseFunction1((uint64_t)(posX + i)) +
                                  noiseFunction2((uint64_t)(posY + j));
            unsigned charIndex = (24-j) * 25 + i;
            charBuffer[charIndex] = (unsigned)(((noiseValue >> 9) * 257) >> (64 - 9));
            noiseValue &= 0x007fffffffffffff;
            fgClrBuffer[charIndex] = (unsigned)((noiseValue & 0x007fffff80000000) >> 31);
            noiseValue &= 0x0000000000ffffff;
            bgClrBuffer[charIndex] = (unsigned)noiseValue;
        }
    }

    cg.WriteCharacters(charBuffer, fgClrBuffer, bgClrBuffer, 0, 0, 25, 25);

    unsigned const* pixelBuffer = cg.GetPixels();
    unsigned cgPixelRows = cg.pixelRows();
    unsigned cgPixelCols = cg.pixelCols();
    for (unsigned j = 0; j < cgPixelRows; j++) {
        for (unsigned i = 0; i < cgPixelCols; i++) {
            //((unsigned*)s->pixels)[(j+2) * s->w + 2 + i] = pixelBuffer[j * cgPixelCols + i];
            ((unsigned*)s->pixels)[(j) * s->w + i] = pixelBuffer[j * cgPixelCols + i];
        }
    }

    SDL_UpdateWindowSurface(w);

    return 0;
}
*/

