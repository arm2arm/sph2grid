/* 
 * File:   Render.cpp
 * Author: arm2arm
 * 
 * Created on January 26, 2010, 5:44 PM
 */

#include "Render.h"
#include <iostream>
#include <CImg.h>
using namespace cimg_library;

CRender::CRender() {
}

CRender::CRender(const CRender& orig) {
}

CRender::~CRender() {
}

void CRender::DoRender(float ***vol3d, int GRID) {

    CImg<float> image(GRID, GRID, 1, 3);

    for (int i = 0; i < GRID - 1; i++)
        for (int j = 0; j < GRID - 1; j++)
            for (int k = 0; k < GRID - 1; k++) {
                float pix = vol3d[i][j][k];
                image(i, j) += pix;
                image(i, j, 0, 1) = image(i, j);
                image(i, j, 2) = image(i, j);

            }



    image = image.normalize(0, 255) + 1;
    (image.log10().normalize(0, 255)).display();



}

