////////////////////////////////////////
//
//  Character graphics for SDL
//
//  2020 October 8
//
////////////////////////////////////////


struct SDLCharGraphicsPane {

    unsigned mNumCharCols;
    unsigned mNumCharRows;
    unsigned mNumCharMarginPixelCols;
    unsigned mNumCharMarginPixelRows;
    unsigned mNumPixelColsPerCharCol;
    unsigned mNumPixelRowsPerCharRow;
    unsigned mCharSizePixelCols;
    unsigned mCharSizePixelRows;
    unsigned mNumPixelCols;
    unsigned mNumPixelRows;
    unsigned mNumFontBitmaps;
    unsigned mPosX;
    unsigned mPosY;
    bool mReticle;

    unsigned* mChars;
    unsigned* mFGColors;
    unsigned* mBGColors;
    unsigned** mFontBitmaps;

    void Clear(unsigned* pixels, unsigned numWindowCols);

    SDLCharGraphicsPane(unsigned charCols, unsigned charRows,
                        unsigned charSizeCols, unsigned charSizeRows,
                        unsigned horizPixelSpace, unsigned vertPixelSpace,
                        unsigned* fontSpace, unsigned fontSpaceCols, unsigned fontSpaceRows,
                        unsigned posX, unsigned posY, bool reticle);
    ~SDLCharGraphicsPane();
    
};


class SDLCharGraphics {

    void UpdatePixelBuffer(SDLCharGraphicsPane* pane);

public:

    unsigned mPixelRows;
    unsigned mPixelCols;
    unsigned* mPixels;

    SDLCharGraphics(unsigned numPixelRows, unsigned numPixelCols);
    virtual ~SDLCharGraphics();

    int WriteCharacters(SDLCharGraphicsPane* pane, unsigned* chars, 
                        unsigned* fgColors, unsigned* bgColors,
                        int x, int y, unsigned w, unsigned h);
    unsigned const* GetPixels();

};
