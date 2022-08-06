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
#include <cstring>
#include <vector>


#define kMainWindowCols     1012
#define kMainWindowRows      512
#define kBitsPerStellarSquare  3
#define kISMapColumns         25
#define kISMapRows            25
#define kIPMapColumns         25
#define kIPMapRows            25
#define kTerrainBuffColumns 4000
#define kTerrainBuffRows    2000
#define kInfoColumns          27
#define kInfoRows             25
#define kMaxPlanetaryLevel     6

// These need to be >= any of the other character Columns and Rows #defines.
// This rule excludes the #defines that refer to number of pixel columsn and rows.
#define kMaxColumns          50
#define kMaxRows             50


static constexpr uint32_t gOrigin = uint32_t(1) << 31;


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
    uint32_t mX;
    uint32_t mY;
    uint32_t mStepSize;
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
static float gTerrainBuffer[kTerrainBuffColumns*kTerrainBuffRows];
static double gWaterline = 0.0;


//GameState* 
void HandleKeyDown(SDL_Event& sdlEvent, SDL_Window* win, SDLCharGraphics& charGraphics);
void ClearCharGraphicsBuffer(unsigned* buffer, unsigned bufferSize);
void FillColorBuffer(unsigned* buffer, unsigned bufferSize, unsigned color);
void CopyStrToIntBuffer(unsigned* buffer, char* str, unsigned maxlen);
unsigned TempToRGB(double temp, double extraSaturation);
unsigned LumToSymbol(double lum);
unsigned PlanetMassToSymbol(double massInKg);
unsigned LifeRGB(double life);
unsigned PlanetSymbolRGB(Planet* planet, StarSystem* starSystem);
//vector<unsigned> PlanetPalette(uint32_t seed, unsigned size);
void FormatSystemInfo(unsigned* buffer, uint32_t posX, uint32_t posY, StarSystem* starSystem);
void TopComponents(vector<double>& v, int& first, int& second, int& third, double& total);
void FormatPlanetComponents(char* line, vector<double>& abundances, Planet* planet);
void FormatPlanetInfo(unsigned* buffer, uint32_t posX, uint32_t posY, Planet* planet);
void PaintTerrain(SDL_Window* w);
int UpdateGraphics(SDLCharGraphics& cg, SDL_Window* win);
inline float PerlinNoise(uint32_t x, uint32_t y, unsigned fracbits);
inline float PerlinNoiseXPeriodic(uint32_t x, uint32_t y, unsigned fracbits, 
                                  uint32_t xOrigin, unsigned periodbits);
inline double PerlinNoise(double x, double y);
inline double PerlinNoise(double x, double y, double z);
inline double Hermite(double x);
inline float Hermite(float x);
//inline uint64_t ToFixedPoint12p52(double x);
inline uint32_t ToFixedPoint9p23(float x);
inline uint64_t NoiseFunction1(uint64_t x);
inline uint64_t NoiseFunction2(uint64_t x);
inline uint32_t NoiseFunction1(uint32_t x);
inline uint32_t NoiseFunction2(uint32_t x);


int main(int argc, char** argv)
{
    //AllocConsole();
    //HANDLE theConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    //freopen("CON", "w", stdout);

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* win = 
        SDL_CreateWindow("Space Game", 40, 40, kMainWindowCols, kMainWindowRows, 0);

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

    SDLCharGraphics charGraphics(kMainWindowRows, kMainWindowCols);

    gInterstellarState->mMapPane = 
        new SDLCharGraphicsPane(kISMapColumns, kISMapRows, 16, 16, 2, 2, 
                                fontSheetPixels, 16, 22, 0, 0, true);
    gInterstellarState->mInfoPane =
        new SDLCharGraphicsPane(kInfoColumns, kInfoRows, 14, 16, 2, 2, 
                                fontSheetPixels2, 16, 22, 512, 0, false);
    gInterstellarState->mStarSystem = nullptr;
    gInterstellarState->mPlanet = nullptr;
    gInterstellarState->mX = 0;
    gInterstellarState->mY = 0;
    gInterstellarState->mStepSize = uint32_t(1) << kBitsPerStellarSquare;
    gInterstellarState->mLevel = 0;
    gInterstellarState->mLastLevel = 0;
    gInterstellarState->mLastState = nullptr;
    gInterstellarState->mID = GameState::kInterstellar;

    gInterplanetaryState->mMapPane = 
        new SDLCharGraphicsPane(kIPMapColumns, kIPMapRows, 16, 16, 2, 2, 
                                fontSheetPixels, 16, 22, 0, 0, true);
    gInterplanetaryState->mInfoPane =
        new SDLCharGraphicsPane(kInfoColumns, kInfoRows, 14, 16, 2, 2, 
                                fontSheetPixels2, 16, 22, 512, 0, false);
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
                                fontSheetPixels2, 16, 22, 512, 0, false);
    gPlanetaryState->mStarSystem = nullptr;
    gPlanetaryState->mPlanet = nullptr;
    gPlanetaryState->mX = 0;
    gPlanetaryState->mY = 0;
    gPlanetaryState->mStepSize = 1 << kMaxPlanetaryLevel;
    gPlanetaryState->mLevel = 0;
    gPlanetaryState->mLastLevel = 0;
    gPlanetaryState->mLastState = nullptr;
    gPlanetaryState->mID = GameState::kPlanetary;
    
    gCurrentState = gInterstellarState;

    gInterstellarState->mX = gOrigin;
    gInterstellarState->mY = gOrigin;
    UpdateGraphics(charGraphics, win);

    SDL_Event sdlEvent;
    do {
        gCurrentState->mLastState = gCurrentState;
        if (0 != SDL_PollEvent(&sdlEvent)) {
            if (SDL_KEYDOWN == sdlEvent.type) {
                HandleKeyDown(sdlEvent, win, charGraphics);
            }
            else if (SDL_WINDOWEVENT == sdlEvent.type &&
                     (SDL_WINDOWEVENT_EXPOSED == sdlEvent.window.event ||
                      SDL_WINDOWEVENT_MAXIMIZED == sdlEvent.window.event ||
                      SDL_WINDOWEVENT_RESTORED == sdlEvent.window.event ||
                      SDL_WINDOWEVENT_FOCUS_GAINED == sdlEvent.window.event)) {
                UpdateGraphics(charGraphics, win);
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


void HandleKeyDown(SDL_Event& sdlEvent, SDL_Window* win, SDLCharGraphics& charGraphics)
{
    uint32_t d = (gCurrentState->mLevel > 0) ?
                 gCurrentState->mStepSize >>  gCurrentState->mLevel :
                 gCurrentState->mStepSize << -gCurrentState->mLevel;
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
            else if (gCurrentState->mLevel < kMaxPlanetaryLevel) {
                gCurrentState->mLevel++;
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
                    gCurrentState->mLevel--;
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
        UpdateGraphics(charGraphics, win);
    }
    gCurrentState = (SDLK_ESCAPE == sdlEvent.key.keysym.sym) ? nullptr : gCurrentState;
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
    else if (planet->mHydCover > 0.25) {
        return 0x003636ff;
    }
    else {
        return 0x00845c18;
    }
}


void FormatSystemInfo(unsigned* buffer, uint32_t posX, uint32_t posY, StarSystem* starSystem)
{
    char line[32];
    sprintf(line, "%010u x %010u", posX, posY);
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


void TopComponents(vector<double>& v, int& first, int& second, int& third, double& total)
{
	first = -1;
	second = -1;
	third = -1;
	total = 0.0;
	double x1 = -1e+300;
	double x2 = -1e+300;
	double x3 = -1e+300;
	for (int i = 0; i < v.size(); ++i) {
		double x = v[i];
		total += x;
		if (x > x3) {
			if (x > x2) {
				x3 = x2;
				third = second;
				if (x > x1) {
					x2 = x1;
					second = first;
					x1 = x;
					first = i;
				}
				else {
					x2 = x;
					second = i;
				}
			}
			else {
				x3 = x;
				third = i;
			}
		}
	}
}


void FormatPlanetComponents(char* line, vector<double>& abundances, Planet* planet)
{
    const auto& comps = planet->mVolComponents; // Just to have a short variable name
    int first, second, third;
    double total, concSum;
    string* form;
    char* cptr;
    TopComponents(abundances, first, second, third, total);
    line[0] = '\0';
    if (total > 0.0) {
        concSum = 0.0;
        form = &kTheComponentFormulae[comps[first]->mComponentIndex];
        cptr = line;
        strcpy(line, form->c_str());
        cptr += form->length();
        concSum += first / total;
        if (concSum < 0.99) {
            form = &kTheComponentFormulae[comps[second]->mComponentIndex];
            sprintf(cptr, " %s", form->c_str());
            cptr += form->length() + 1;
            concSum += second / total;
            if (concSum < 0.99) {
                sprintf(cptr, " %s", 
                        kTheComponentFormulae[comps[third]->mComponentIndex].c_str());
            }
        }
    }
}


void FormatPlanetInfo(unsigned* buffer, uint32_t posX, uint32_t posY, Planet* planet)
{
    char line[32];
    sprintf(line, "%010u x %010u", posX, posY);
    CopyStrToIntBuffer(buffer, line, kInfoColumns);
    
    if (planet != nullptr) {
		sprintf(line, "Temp %.0lf K", planet->mTemperature);
		CopyStrToIntBuffer(buffer+kInfoColumns, line, kInfoColumns);
		sprintf(line, "Pres %.3g atm", planet->mPressure);
		CopyStrToIntBuffer(buffer+kInfoColumns*2, line, kInfoColumns);
		
		const auto& comps = planet->mVolComponents; // Just to have a short variable name
		vector<double> abundances(comps.size());
		
		for (int i = 0; i < comps.size(); ++i) {
			abundances[i] = comps[i]->mAtmAbundance;
		}
        FormatPlanetComponents(line, abundances, planet);
        CopyStrToIntBuffer(buffer+kInfoColumns*3, line, kInfoColumns);
		
        strcpy(line, "Seas: ");
		for (int i = 0; i < comps.size(); ++i) {
			abundances[i] = comps[i]->mHydAbundance;
		}
        FormatPlanetComponents(line+6, abundances, planet);
        CopyStrToIntBuffer(buffer+kInfoColumns*4, line, kInfoColumns);
		
        strcpy(line, "Clouds: ");
		for (int i = 0; i < comps.size(); ++i) {
			abundances[i] = comps[i]->mCldAbundance;
		}
        FormatPlanetComponents(line+8, abundances, planet);
        CopyStrToIntBuffer(buffer+kInfoColumns*5, line, kInfoColumns);
		
        strcpy(line, "Frost: ");
		for (int i = 0; i < comps.size(); ++i) {
			abundances[i] = comps[i]->mFstAbundance;
		}
        FormatPlanetComponents(line+7, abundances, planet);
        CopyStrToIntBuffer(buffer+kInfoColumns*6, line, kInfoColumns);
	}
}


// Gives a value, the partial power series sum of which equals 1
// i.e. returns x such that x + x/2 + ... + x/2^(n-1) = 1
inline float UnityFirstTerm(unsigned numTerms)
{
    return float(1 << (numTerms-1)) / float((1 << numTerms) - 1);
}


void PaintTerrain(SDL_Window* w)
{
    vector<float> analysisBuffer;
    vector<uint32_t> warpXBuffer;
    vector<uint32_t> warpYBuffer;
    vector<float> multiBuffer;
    vector<float> powerBuffer;
    const uint32_t terrainColsBits = 11;
    const uint32_t terrainScaleBits = terrainColsBits - 2;
    uint32_t terrainCols = 1 << terrainColsBits;
    uint32_t terrainRows = 1 << (terrainColsBits-1);
    float terrainScale = float(1 << terrainScaleBits);
    if (gPlanetaryState->mLastState != gPlanetaryState) {
        uint32_t paramSeed = 
            NoiseFunction1(uint32_t(gCurrentPlanet->mSemimajorAxis * float(uint64_t(1)<<32)));
        uint32_t planetx = (paramSeed >> terrainScaleBits) << terrainScaleBits;
        paramSeed = NoiseFunction1(paramSeed);
        uint32_t planety = paramSeed;
        paramSeed = NoiseFunction1(paramSeed);
        uint32_t warpxx = (paramSeed >> terrainScaleBits) << terrainScaleBits;
        paramSeed = NoiseFunction1(paramSeed);
        uint32_t warpxy = paramSeed;
        paramSeed = NoiseFunction1(paramSeed);
        uint32_t warpyx = (paramSeed >> terrainScaleBits) << terrainScaleBits;
        paramSeed = NoiseFunction1(paramSeed);
        uint32_t warpyy = paramSeed;
        paramSeed = NoiseFunction1(paramSeed);
        uint32_t multix = (paramSeed >> terrainScaleBits) << terrainScaleBits;
        paramSeed = NoiseFunction1(paramSeed);
        uint32_t multiy = paramSeed;
        float warpStrength = float(paramSeed) / float(uint64_t(1)<<32) * 1.2f + 0.15f;
        float scalePerOctave = 0.5f;
        unsigned numOctaves = 1;
        float initOctaveScale = warpStrength * UnityFirstTerm(numOctaves);
        for (uint32_t j = 0; j < terrainRows; ++j) {
            for (uint32_t i = 0; i < terrainCols; ++i) {
                uint32_t xbufval = 0;
                uint32_t ybufval = 0;
                float octaveScale = initOctaveScale;
                for (uint32_t n = 0; n < numOctaves; ++n) {
                    float wx = PerlinNoiseXPeriodic((i<<n)+warpxx, (j<<n)+warpxy, 
                                                    terrainScaleBits, 
                                                    warpxx, terrainColsBits+n);
                    float wy = PerlinNoiseXPeriodic((i<<n)+warpyx, (j<<n)+warpyy, 
                                                    terrainScaleBits, 
                                                    warpyx, terrainColsBits+n);
                    wx *= wx;
                    wy *= wy;  // Yes, this transformation
                    wx -= 0.5; // is asymmetric around 0
                    wy -= 0.5;
                    float totalScale = 2.0f * octaveScale * terrainScale * warpStrength;
                    xbufval += uint32_t(wx * totalScale + 0.5);
                    ybufval += uint32_t(wy * totalScale + 0.5);
                    octaveScale *= scalePerOctave;
                }
                warpXBuffer.push_back(xbufval);
                warpYBuffer.push_back(ybufval);
            }
        }
        scalePerOctave = 0.5f;
        numOctaves = 1;
        initOctaveScale = UnityFirstTerm(numOctaves);
        for (uint32_t j = 0; j < terrainRows; ++j) {
            for (uint32_t i = 0; i < terrainCols; ++i) {
                unsigned index = j*terrainCols + i;
                float h = 0.0f;
                float octaveScale = initOctaveScale;
                for (uint32_t n = 0; n < numOctaves; ++n) {
                    h += PerlinNoiseXPeriodic(((i<<n)+warpXBuffer[index]) + planetx, 
                                              ((j<<n)+warpYBuffer[index]) + planety, 
                                              terrainScaleBits, planetx, terrainColsBits + n) *
                         octaveScale;
                    octaveScale *= scalePerOctave;
                }
                // This modulates the scalePerOctave of the terrain. Range 0.35-0.7
                multiBuffer.push_back(h*0.35f + 0.35f);
            }
        }
        for (uint32_t j = 0; j < terrainRows; ++j) {
            for (uint32_t i = 0; i < terrainCols; ++i) {
                unsigned index = j*terrainCols + i;
                float h = 0.0f;
                float octaveScale = initOctaveScale;
                for (uint32_t n = 0; n < numOctaves; ++n) {
                    h += PerlinNoiseXPeriodic(((i<<n)+warpXBuffer[index]) + planetx, 
                                              ((j<<n)+warpYBuffer[index]) + planety, 
                                              terrainScaleBits, planetx, terrainColsBits + n) *
                         octaveScale;
                    octaveScale *= scalePerOctave;
                }
                powerBuffer.push_back(h);
            }
        }
        numOctaves = terrainColsBits - 1;
        // UnityFirstTerm is tuned for powers of 1/2, but somehow it's perfect for
        // our purposes here also, where scale per buffer is not generally 1/2. I
        // haven't bothered to figure out why this is.
        initOctaveScale = UnityFirstTerm(numOctaves);
        for (uint32_t j = 0; j < terrainRows; ++j) {
            for (uint32_t i = 0; i < terrainCols; ++i) {
                unsigned index = j*terrainCols + i;
                float h = 0.0f;
                float octaveScale = initOctaveScale;
                scalePerOctave = multiBuffer[index];
                float power = powerBuffer[index];
                for (uint32_t n = 0; n < numOctaves; ++n) {
                    float oct = PerlinNoiseXPeriodic((i<<n) + warpXBuffer[index] + planetx, 
                                                     (j<<n) + warpYBuffer[index] + planety, 
                                                     terrainScaleBits, planetx, 
                                                     terrainColsBits + n);
                    float octcubed = oct * oct * oct;
                    oct = (1.0f-power) * oct + power * octcubed;
                    oct *= octaveScale;
                    h += oct;
                    octaveScale *= scalePerOctave;
                }
                gTerrainBuffer[index] = h;
                analysisBuffer.push_back(h);
            }
        }
        sort(analysisBuffer.begin(), analysisBuffer.end());
        gWaterline = analysisBuffer[unsigned(floor((analysisBuffer.size()-1) *
                                                   sqrt(gCurrentPlanet->mHydCover)))];
    }

    SDL_Surface* surf = SDL_GetWindowSurface(w);
    unsigned* pixelsPtr = ((unsigned*)surf->pixels);

    int maxIndex = 1 << (2*terrainColsBits - 1);
    for (int j = 0; j < 512; ++j) {
        for (int i = 0; i < 512; ++i) {
            // j: [0, 128, 256, 384, 512] ==> 
            // [rows*3/2, rows, rows/2, 0, -rows/2]  {level 0}
            // [rows, rows*3/4, rows/2, rows/4, 0]   {level 1}
            // [rows*3/4, rows*5/8, rows/2, rows*3/8, rows/4]  {level 2}
            // rows/2 + [rows>>level, rows/2>>level, 0, -rows/2>>level, -rows>>level] =
            // rows/2 + rows>>level + 
            //  [0, -rows/2>>level, -rows>>level, -rows*3/2>>level, -rows*2>>level] =
            // rows/2 + rows>>level - ((rows*j)>>8)>>level
            const int lev = gCurrentState->mLevel; // Just for a short variable name
            int index = ((((i << (terrainColsBits-9)) >> lev) -
                          (terrainCols/2 >> lev) + gCurrentState->mX) & terrainCols-1) + 
                        (((terrainRows>>1) + (terrainRows>>lev) - (terrainRows*j >> (lev+8)) + 
                          gCurrentState->mY) << terrainColsBits);
            if (index < 0 || index >= maxIndex) {
                pixelsPtr[j * surf->w + i] = 0;
            }
            else if (gTerrainBuffer[index] > gWaterline) {
                unsigned s = gTerrainBuffer[index] * 256.0;
                s = (s < 0x80000000) ? s : 0;
                s = (s <= 255) ? s : 255;
                pixelsPtr[j * surf->w + i] = 0x00010000 * s;
                pixelsPtr[j * surf->w + i] += 0x00000100 * (s>>1);
                pixelsPtr[j * surf->w + i] += 0x00000001 * (s>>2);
            }
            else {
                pixelsPtr[j * surf->w + i] = 0x00004080;
            }
        }
    }
}


int UpdateGraphics(SDLCharGraphics& cg, SDL_Window* w)
{
    SDL_Surface* surf = SDL_GetWindowSurface(w);
    
    if (GameState::kInterstellar == gCurrentState->mID) {
        gSystemUnderReticle = nullptr;
        unsigned halfstep = gInterstellarState->mStepSize >> 1;
        for (int j = -12; j < 13; j++) {
            for (int i = -12; i < 13; i++) {
                // Adding half a step so the value isn't exactly 0.5 at integer x and y.
                float noiseValue = 
                    PerlinNoise(gCurrentState->mX + (unsigned(i) << kBitsPerStellarSquare) +
                                halfstep, 
                                gCurrentState->mY + (unsigned(j) << kBitsPerStellarSquare) +
                                halfstep, 6);
                float maxNoiseValue = 1.0f - 1.0f/float(1 << 24);
                noiseValue = (noiseValue < 0.0) ? 0.0 :
                             ((noiseValue > maxNoiseValue) ? maxNoiseValue : noiseValue);
                float starNoise = noiseValue * (1 << 11);
                float starFrac = starNoise - floor(starNoise);
                float starFrequency = (63.0f * pow(noiseValue, 8.0f) + 1.0f) / 192.0f;
                int index = (12-j)*25 + i + 12;
                cgBuffer[index] = 32;
                fgClrBuffer[index] = 0xffffff;
                bgClrBuffer[index] = unsigned(noiseValue * 256.0f) | (32 << 16);
                StarSystem* tempSystem = nullptr;
                bool center = ((i == 0) && (j == 0));
                if (starFrac <= starFrequency) {
                    uint32_t fpNoiseValue = ToFixedPoint9p23(noiseValue);
                    tempSystem = new StarSystem(fpNoiseValue, !center);
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
        FormatSystemInfo(cgBuffer, gCurrentState->mX, gCurrentState->mY, gCurrentSystem);
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
        double orbitScale = 1.0;
        if (numPlanets > 1) {
            orbitScale = max(1.01 / log(minRatio), 90.0 / log(maxSA/minSA));
        }
        vector<uint32_t> planetYs;
        uint32_t minY = -1; // It'll underflow to the max unsigned
        for (unsigned i = 0; i < numPlanets; ++i) {
            Planet* p = gCurrentSystem->mPlanets[i];
            uint32_t y = uint32_t(gOrigin - 5) - 
                         uint32_t(floor(orbitScale * (log(p->mSemimajorAxis)-log(minSA))));
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
        uint32_t screenX = -gCurrentState->mX + (gOrigin + 12);
        uint32_t screenY =  gCurrentState->mY - (gOrigin - 12);
        if (screenX < 25 && screenY < 25) {
            unsigned index = screenX + 25 * screenY;
            cgBuffer[index] = 42;
            fgClrBuffer[index] = TempToRGB(gCurrentSystem->mStar->mTemperature, 0.5);
            bgClrBuffer[index] = 0x00000000;
        }
        // Planets
        int planetUnderReticle = -1;
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
                    planetUnderReticle = i;
                }
            }
        }
        cg.WriteCharacters(gCurrentState->mMapPane, cgBuffer, fgClrBuffer, bgClrBuffer, 0, 0, 
                           kIPMapColumns, kIPMapRows);
        FillColorBuffer(bgClrBuffer, kMaxColumns*kMaxRows, 0);
        if (planetUnderReticle != -1) {
            for (int i = 0; i < kInfoColumns; ++i) {
                bgClrBuffer[(planetUnderReticle+2)*kInfoColumns + i] = 0x05005f;
            }
        }
        cg.WriteCharacters(gInterstellarState->mInfoPane, nullptr, nullptr, bgClrBuffer, 
                           0, 0, kInfoColumns, kInfoRows);
    }
    else {
        PaintTerrain(w); // This paints directly to the screen as a side effect
        ClearCharGraphicsBuffer(cgBuffer, kMaxColumns*kMaxRows);
        FillColorBuffer(fgClrBuffer, kMaxColumns*kMaxRows, 0xffffff);
        FillColorBuffer(bgClrBuffer, kMaxColumns*kMaxRows, 0);
        FormatPlanetInfo(cgBuffer, gCurrentState->mX, gCurrentState->mY, gCurrentPlanet);
        cg.WriteCharacters(gCurrentState->mInfoPane, cgBuffer, fgClrBuffer, bgClrBuffer, 0, 0, 
                           kInfoColumns, kInfoRows);
    }

    unsigned paintColsStart = (GameState::kPlanetary == gCurrentState->mID) ?
                              gPlanetaryState->mInfoPane->mPosX : 0;
    unsigned const* pixelBuffer = cg.GetPixels();
    for (unsigned j = 0; j < kMainWindowRows; j++) {
        for (unsigned i = paintColsStart; i < kMainWindowCols; i++) {
            ((unsigned*)surf->pixels)[j * surf->w + i] = 
                pixelBuffer[j * kMainWindowCols + i];
        }
    }

    SDL_UpdateWindowSurface(w);

    return 0;
}


float PerlinNoise(uint32_t x, uint32_t y, unsigned fracbits)
{
    uint32_t xint = x >> fracbits;
    uint32_t yint = y >> fracbits;
    uint32_t xp1 = (x + (1 << fracbits)) >> fracbits;
    uint32_t yp1 = (y + (1 << fracbits)) >> fracbits;
    uint32_t ymod = yint * 1103515245 + 12345;
    uint32_t yp1mod = yp1 * 1103515245 + 12345;
    uint32_t ll = NoiseFunction1(xint ^ ymod);
    uint32_t lr = NoiseFunction1( xp1 ^ ymod);
    uint32_t ul = NoiseFunction1(xint ^ yp1mod);
    uint32_t ur = NoiseFunction1( xp1 ^ yp1mod);
    
    // 23 bits in single precision floating point mantissa
    // No point in keeping more
    const float phaseNormalizer = 3.14159267f * 2.0f / 8388608.0f;
    float llphase = float(ll & 0x007fffff) * phaseNormalizer;
    float lrphase = float(lr & 0x007fffff) * phaseNormalizer;
    float ulphase = float(ul & 0x007fffff) * phaseNormalizer;
    float urphase = float(ur & 0x007fffff) * phaseNormalizer;

    unsigned fracmask = (1 << fracbits) - 1;    
    float fracNormalizer = 1.0f / float(fracmask);
    float xfrac = float(x & fracmask) * fracNormalizer;
    float yfrac = float(y & fracmask) * fracNormalizer;

    float lldot =  xfrac      *cos(llphase) +  yfrac      *sin(llphase);
    float lrdot = (xfrac-1.0f)*cos(lrphase) +  yfrac      *sin(lrphase);
    float uldot =  xfrac      *cos(ulphase) + (yfrac-1.0f)*sin(ulphase);
    float urdot = (xfrac-1.0f)*cos(urphase) + (yfrac-1.0f)*sin(urphase);
    
    float hermitex = Hermite(xfrac);
    float hermitey = Hermite(yfrac);

    float linterp = (1.0f-hermitex)*lldot + hermitex*lrdot;
    float uinterp = (1.0f-hermitex)*uldot + hermitex*urdot;
    float interp = (1.0f-hermitey)*linterp + hermitey*uinterp;

    float result = interp * 0.70710678f + 0.5f;
    return result;
}


#define kArccosPiOver8 0.9238795f
#define kArcsinPiOver8 0.3826834f
inline float GradDotProduct(unsigned n, float x, float y)
{
    switch (n)
    {
        case 0:
        return kArccosPiOver8*x + kArcsinPiOver8*y;
        break;
        
        case 1:
        return kArccosPiOver8*y + kArcsinPiOver8*x;
        break;
        
        case 2:
        return kArccosPiOver8*x - kArcsinPiOver8*y;
        break;
        
        case 3:
        return kArccosPiOver8*y - kArcsinPiOver8*x;
        break;
      
        case 4:
        return -kArccosPiOver8*x + kArcsinPiOver8*y;
        break;
        
        case 5:
        return -kArccosPiOver8*y + kArcsinPiOver8*x;
        break;
        
        case 6:
        return -kArccosPiOver8*x - kArcsinPiOver8*y;
        break;
        
        case 7:
        return -kArccosPiOver8*y - kArcsinPiOver8*x;
        break;
                
        default:
        return 0;
        break;
    }
}


inline
float PerlinNoiseXPeriodic(uint32_t x, uint32_t y, unsigned fracbits, 
                           uint32_t xOrigin, unsigned periodbits)
{
    // Throw away last bits of xOrigin because it must be an integer
    // with respect to this function's sense of fixed point arithmetic
    xOrigin >>= fracbits;
    xOrigin <<= fracbits;
    
    uint32_t periodMask = (1 << periodbits) - 1;
     
    x -= xOrigin;
    uint32_t xp1 = x + (1 << fracbits);
    x &= periodMask;
    xp1 &= periodMask;
    x += xOrigin;
    xp1 += xOrigin;

    uint32_t xint = x >> fracbits;
    xp1 >>= fracbits;
    uint32_t yint = y >> fracbits;
    uint32_t ymod = yint * 1103515245 + 12345;
    uint32_t ll = NoiseFunction1(xint ^  ymod);
    uint32_t lr = NoiseFunction1(xp1  ^  ymod);
    uint32_t ul = NoiseFunction1(xint ^ (ymod+1103515245));
    uint32_t ur = NoiseFunction1(xp1  ^ (ymod+1103515245));
    
    uint32_t fracMask = (1 << fracbits) - 1;
    float fracNormalizer = 1.0f / float(1 << fracbits);
    float xfrac = float(x & fracMask) * fracNormalizer;
    float yfrac = float(y & fracMask) * fracNormalizer;

    // gi stands for gradient index
    unsigned llgi = ll >> 29;
    unsigned lrgi = lr >> 29;
    unsigned ulgi = ul >> 29;
    unsigned urgi = ur >> 29;

    float xfracm1 = xfrac - 1.0;
    float yfracm1 = yfrac - 1.0;

    float lldot = GradDotProduct(llgi, xfrac,   yfrac);
    float lrdot = GradDotProduct(lrgi, xfracm1, yfrac);
    float uldot = GradDotProduct(ulgi, xfrac,   yfracm1);
    float urdot = GradDotProduct(urgi, xfracm1, yfracm1);
    
    float hermitex = Hermite(xfrac);
    float hermitey = Hermite(yfrac);

    float linterp = (1.0f-hermitex)*lldot + hermitex*lrdot;
    float uinterp = (1.0f-hermitex)*uldot + hermitex*urdot;
    float interp = (1.0f-hermitey)*linterp + hermitey*uinterp;

    float result = interp * 0.70710678f + 0.5f;
    return result;
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


inline uint32_t NoiseFunction1(uint32_t x)
{
    uint32_t result = x * 1103515245 + 12345;
    result ^= result << 13;
    result ^= result >> 7;
    result ^= result << 17;
    return result;
}


inline uint32_t NoiseFunction2(uint32_t x)
{
    uint32_t result = x * 6680178296815197433 + 7046029254386353087;
    result ^= result << 13;
    result ^= result >> 7;
    result ^= result << 17;
    return result;
}

// x should be in [0,1]
inline double Hermite(double x)
{
    double xsquared = x * x;
    return x * xsquared * (6.0*xsquared - 15.0*x + 10.0);
}

// x should be in [0,1]
inline float Hermite(float x)
{
    float xsquared = x * x;
    return x * xsquared * (6.0f*xsquared - 15.0f*x + 10.0f);
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

// x should be in [0,1]
inline uint32_t ToFixedPoint9p23(float x)
{
    uint32_t crosscast = *((uint32_t*)(&x));
    uint32_t mantissa = 0x00800000 | (crosscast & 0x007fffff);
    int32_t exponent = int32_t((0x7f800000 & crosscast) >> 23) - 127;
    exponent = (exponent > 0) ? 0 : exponent;
    exponent = (exponent < -23) ? -23 : exponent;
    uint32_t result = mantissa >> -exponent;
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

