#define SDL_MAIN_HANDLED                                                // windows stuff

#include "../include/params.h"

#include "../include/Dirac/dirac.h"

#include "../include/Griddy/griddy.h"
#include "../include/Griddy/SDL2/SDL.h"
#include "../include/Griddy/SDL2/SDL_image.h"
#include <iostream>
#include <cmath>

int main(int argc, char* args[])
{
    Base       xBase = Base(START_X, SIZE_X, END_X);
    Base       yBase = Base(START_Y, SIZE_Y, END_Y);
    Basis2*    basis = new Basis2(xBase, yBase)    ;
    Scalar2*   u     = new Scalar2(basis)          ;
    WaveFunc2* psi   = new WaveFunc2(basis)        ;

    float          potnVals[SIZE_X][SIZE_Y];
    Complex        probAmps[SIZE_X][SIZE_Y];
    griddy::Vertex vertices[SIZE_X][SIZE_Y];
    
    griddy::Window window = griddy::Window(
        WINDOW_TITLE ,
        WINDOW_WIDTH ,
        WINDOW_HEIGHT
    );

    griddy::Grid grid = griddy::Grid(&window);

    SDL_Rect rectGrid = SDL_Rect();
    rectGrid.x = 0                ;
    rectGrid.y = 0                ;
    rectGrid.w = WINDOW_WIDTH     ;
    rectGrid.h = WINDOW_HEIGHT    ;

    grid.setVertices(&vertices[0][0]);
    grid.setPosition(rectGrid);
    grid.setSize(SIZE_X, SIZE_Y);

    u  ->setValues(&potnVals[0][0]);
    psi->setValues(&probAmps[0][0]);

    psi->setNorm(MAX_COLOR);
    psi->setMass(1.0f);

    for (int i = 0; i < SIZE_X; i++)
        for (int j = 0; j < SIZE_Y; j++)
        {
            float x = xBase(i) - INIT_X;
            float y = yBase(j) - INIT_Y;

            probAmps[i][j] = Real(exp((-x*x-y*y) / DEV)) * cis(MOMENTUM_X * x + MOMENTUM_Y * y);
            potnVals[i][j] = 0.5f;
            vertices[i][j] = griddy::Vertex(i, j);
        }

    SDL_Event event;
    while (window.isRunning)
    {
        while ( SDL_PollEvent(&event) )
            if (event.type == SDL_QUIT)
                window.destroy();

        window.clear();

        psi->evolve(DT, u);
        float factor = psi->prbFactor(); 
        for (int i = 0; i < SIZE_X; i++)
            for (int j = 0; j < SIZE_Y; j++)
            {
                float thisProb = factor * psi->prob(i, j, false);
                griddy::color thisColor = griddy::Color(0x00, thisProb, 0x00);
                vertices[i][j].setColor(thisColor);
            }
        grid.render();

        window.display();
    }

    SDL_Quit();                                                         // goodbye
    return 0;
}