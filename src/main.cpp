#define SDL_MAIN_HANDLED                                                // windows stuff

#include "../include/params.h"

#include "../include/Dirac/dirac.h"

#include "../include/Griddy/griddy.h"
#include "../include/Griddy/SDL2/SDL.h"
#include <iostream>
#include <cmath>

int main(int argc, char* args[])
{
    griddy::testVideo();                                    // test SDL

    Base*      xBase = new Base(START_X, SIZE_X, END_X);    // x base
    Base*      yBase = new Base(START_Y, SIZE_Y, END_Y);    // y base
    Basis2*    basis = new Basis2(xBase, yBase)        ;    // 2D basis
    Scalar2*   field = new Scalar2(basis)              ;    // potential field AKA the slits
    WaveFunc2* psi   = new WaveFunc2(basis)            ;    // the wave function AKE the "particle"

    float         fieldVals[SIZE_X][SIZE_Y];                // values for the potential field
    Complex        probAmps[SIZE_X][SIZE_Y];                // values for the wave function
    griddy::Vertex vertices[SIZE_X][SIZE_Y];                // the "pixels" of the screen
    
    griddy::Window window = griddy::Window(
        WINDOW_TITLE ,
        WINDOW_WIDTH ,
        WINDOW_HEIGHT
    );

    griddy::Grid grid = griddy::Grid(&window);

    SDL_Rect rectGrid = SDL_Rect();                         // position of the grid on the window (null position)
    rectGrid.x = 0                ;
    rectGrid.y = 0                ;
    rectGrid.w = WINDOW_WIDTH     ;
    rectGrid.h = WINDOW_HEIGHT    ;

    grid.setVertices(&vertices[0][0]);
    grid.setPosition(rectGrid);
    grid.setSize(SIZE_X, SIZE_Y);

    field->setValues(&fieldVals[0][0]);
    psi  ->setValues(& probAmps[0][0]);
    psi  ->setNorm  (MAX_COLOR);
    psi  ->setMass  (MASS);

    for (int i = 0; i < SIZE_X; i++)
        for (int j = 0; j < SIZE_Y; j++)
        {
            float x = xBase->coord(i) - INIT_X;
            float y = yBase->coord(j) - INIT_Y;

            probAmps[i][j] = Real(exp((-x*x-y*y) / DEV)) * cis(- MOMENTUM_X * x - MOMENTUM_Y * y);
            vertices[i][j] = griddy::Vertex(i, j);

            if (
                x > SLIT_X &&
                x < SLIT_X + WIDTH_SLIT
            )    fieldVals[i][j] = POTENTIAL;
            else fieldVals[i][j] = 0.0f;
        }

    SDL_Event event;
    while (window.isRunning)
    {
        while ( SDL_PollEvent(&event) )
            if (event.type == SDL_QUIT)
                window.destroy();

        window.clear();

        psi->evolve(DT, field);
        float factor = psi->prbFactor(); 
        for (int i = 0; i < SIZE_X; i++)
            for (int j = 0; j < SIZE_Y; j++)
            {
                if (fieldVals[i][j] == 0.0f)
                {
                    float thisProb = factor * psi->prob(i, j, false);
                    griddy::color thisColor = griddy::Color(0x00, thisProb, 0x00);
                    vertices[i][j].setColor(thisColor);
                }
                else
                {
                    griddy::color thisColor = griddy::Color(0x0a, 0x0a, 0x0a);
                    vertices[i][j].setColor(thisColor);
                }
            }
        grid.render();

        window.display();
    }

    SDL_Quit();                                                         // goodbye
    return 0;
}