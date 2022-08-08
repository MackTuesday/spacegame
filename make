#!/bin/zsh

compile () {
    if [[ $1.cpp -nt $1.o || ! -a $1.o ]]
    then
        echo "Compiling $1.cpp"
        g++ $2 $3 $4 $5 $6 $7 $8 $9 -c $1.cpp -o $1.o
    fi
}

compile lodepng -std=c++14
compile SDLCharGraphics -std=c++14 -g -I/Library/Frameworks/SDL2.framework/Headers
compile game -std=c++14 -O3 -I/Library/Frameworks/SDL2.framework/Headers -I../qrplanets
g++ -o game /Library/Frameworks/SDL2.framework/SDL2 ../qrplanets/qrplanets.a SDLCharGraphics.o lodepng.o game.o
