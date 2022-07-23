#include <cstring>
#include "SDLCharGraphics.h"


SDLCharGraphicsPane::SDLCharGraphicsPane(unsigned charCols, unsigned charRows,
                                         unsigned charSizeCols, unsigned charSizeRows,
                                         unsigned horizPixelMargin, unsigned vertPixelMargin,
                                         unsigned* fontSheet, unsigned fontSheetCols,
                                         unsigned fontSheetRows,
                                         unsigned posX, unsigned posY, bool reticle) :
    mNumCharCols(charCols), mNumCharRows(charRows),
    mNumCharMarginPixelCols(horizPixelMargin), mNumCharMarginPixelRows(vertPixelMargin),
    mCharSizePixelCols(charSizeCols), mCharSizePixelRows(charSizeRows),
    mPosX(posX), mPosY(posY), mReticle(reticle)
{
    mNumPixelColsPerCharCol = mCharSizePixelCols + 2*mNumCharMarginPixelCols;
    mNumPixelRowsPerCharRow = mCharSizePixelRows + 2*mNumCharMarginPixelRows;
    mNumPixelCols = mNumPixelColsPerCharCol * mNumCharCols;// - charSpacePixelCols;
    mNumPixelRows = mNumPixelRowsPerCharRow * mNumCharRows;// - charSpacePixelRows;
    mChars = new unsigned[mNumCharCols * mNumCharRows];
    mFGColors = new unsigned[mNumCharCols * mNumCharRows];
    mBGColors = new unsigned[mNumCharCols * mNumCharRows];

    mNumFontBitmaps = fontSheetCols * fontSheetRows;
    mFontBitmaps = new unsigned*[mNumFontBitmaps];
    for (unsigned i = 0; i < mNumFontBitmaps; i++) {
        mFontBitmaps[i] = new unsigned[mCharSizePixelCols * mCharSizePixelRows];
    }

    for (unsigned j = 0; j < fontSheetRows; j++) {
        for (unsigned i = 0; i < fontSheetCols; i++) {
            for (unsigned y = 0; y < mCharSizePixelRows; y++) {
                for (unsigned x = 0; x < mCharSizePixelCols; x++) {
                    (mFontBitmaps[j * fontSheetCols + i])[y * mCharSizePixelCols + x] =
                        fontSheet[(j * mCharSizePixelRows + y) *
                                  fontSheetCols * mCharSizePixelCols +
                                  i * mCharSizePixelCols + x];
                }
            }
        }
    }
}


SDLCharGraphicsPane::~SDLCharGraphicsPane()
{
    delete[] mChars;
    delete[] mFGColors;
    delete[] mBGColors;
    for (unsigned i = 0; i < mNumFontBitmaps; i++) {
        delete[] mFontBitmaps[i];
    }
    delete[] mFontBitmaps;
}


void SDLCharGraphicsPane::Clear(unsigned* pixels, unsigned numWindowCols)
{
    for (unsigned j = 0; j < mNumPixelRows; j++) {
        for (unsigned i = 0; i < mNumPixelCols; i++) {
            pixels[(j+mPosY) * numWindowCols + i+mPosX] = 0;
        }
    }
}


SDLCharGraphics::SDLCharGraphics(unsigned numPixelRows, unsigned numPixelCols) :
    mPixelRows(numPixelRows), mPixelCols(numPixelCols)
{
    mPixels = new unsigned[numPixelRows*numPixelCols];
}


SDLCharGraphics::~SDLCharGraphics()
{
    delete[] mPixels;
}


// charsCols and charsRows refer to the logical dimensions of the chars array.
// This function writes all of charBuffer if possible. You can't make it write
//  an arbitrary subset of charBuffer.
int SDLCharGraphics::WriteCharacters(SDLCharGraphicsPane* pane, unsigned* charBuffer, 
                                     unsigned* fgClrBuffer, unsigned* bgClrBuffer, 
                                     int destX, int destY, unsigned cbCols, unsigned cbRows)
{
    // Yuck tons of branching. This can be done faster but let's not optimize yet.
    for (int fromY = 0; fromY < (int)cbRows; fromY++) {
        int toY = fromY + destY;
        if (toY < 0 || toY >= pane->mNumCharRows) {
            continue;
        }
        int fromRowStartIndex = fromY * cbCols;
        int toRowStartIndex = toY * pane->mNumCharCols;
        for (int fromX = 0; fromX < (int)cbCols; fromX++) {
            int toU = fromX + destX;
            if (toU < 0 || toU >= pane->mNumCharCols) {
                continue;
            }
            int fromIndex = fromRowStartIndex + fromX;
            int toIndex = toRowStartIndex + toU;
            if (nullptr != charBuffer) {
                pane->mChars[toIndex] = charBuffer[fromIndex];
            }
            if (nullptr != fgClrBuffer) {
                pane->mFGColors[toIndex] = fgClrBuffer[fromIndex];
            }
            if (nullptr != bgClrBuffer) {
                pane->mBGColors[toIndex] = bgClrBuffer[fromIndex];
            }
        }
    }

    UpdatePixelBuffer(pane);

    return 0;
}


void SDLCharGraphics::UpdatePixelBuffer(SDLCharGraphicsPane* pane)
{
    pane->Clear(mPixels, mPixelCols);

    for (unsigned j = 0; j < pane->mNumCharRows; j++) {
        unsigned charsRowStartIndex = j * pane->mNumCharCols;
        unsigned pixelsStartRow = j * pane->mNumPixelRowsPerCharRow + pane->mPosY;
        for (unsigned i = 0; i < pane->mNumCharCols; i++) {
            unsigned charsIndex = charsRowStartIndex + i;
            unsigned pixelsStartCol = i * pane->mNumPixelColsPerCharCol + pane->mPosX;
            for (unsigned y = 0; y < pane->mNumPixelRowsPerCharRow; y++) {
                unsigned bitmapRowStartIndex =
                    (y - pane->mNumCharMarginPixelRows) * pane->mCharSizePixelCols;
                for (unsigned x = 0; x < pane->mNumPixelColsPerCharCol; x++) {
                    unsigned* pixelsPtr = mPixels + (pixelsStartRow+y) * mPixelCols +
                                          pixelsStartCol+x;
                    if (pane->mReticle &&
                        i == pane->mNumCharCols/2 && j == pane->mNumCharRows/2 && 
                        (x == 0 || x == pane->mNumPixelColsPerCharCol-1) &&
                        (y == 0 || y == pane->mNumPixelRowsPerCharRow-1)) {
                        *pixelsPtr = 0x00ffffff;
                    }
                    else 
                    if (x < pane->mNumPixelColsPerCharCol - pane->mNumCharMarginPixelCols &&
                        y < pane->mNumPixelRowsPerCharRow - pane->mNumCharMarginPixelRows &&
                        x >= pane->mNumCharMarginPixelCols && 
                        y >= pane->mNumCharMarginPixelRows) {
                        unsigned bitmapIndex = bitmapRowStartIndex + x -
                                               pane->mNumCharMarginPixelCols;
                        unsigned fontPixel =
                            (pane->mFontBitmaps[pane->mChars[charsIndex]])[bitmapIndex] >> 8;
                        if (fontPixel == 0) {
                            *pixelsPtr = pane->mBGColors[charsIndex];
                        } else {
                            *pixelsPtr = pane->mFGColors[charsIndex];
                        }
                    } 
                    else {
                        *pixelsPtr = pane->mBGColors[charsIndex];
                    }
                }
            }
        }
    }
}


unsigned const* SDLCharGraphics::GetPixels()
{
    return mPixels;
}
