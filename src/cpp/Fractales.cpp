#include <iostream>
#include <cmath>

#include "Fractales.hpp"

FractaleImage::FractaleImage(int m, int n, DomaineMaths domain, int N) : ImageFonctionelSelectionMOOs(m,n,domain), N(N) {
    //Nothing to init
}

void FractaleImage::incN(){
    ++N;

    refreshAll(getCurrentDomaine());
}

void FractaleImage::onDomaineChangePerformed(const DomaineMaths& domainNew){
    N = 10;

    //Repaint everything
    refreshAll(domainNew);
}

FractaleGLImage::FractaleGLImage(FractaleImage* image) : GLImageFonctionelSelections(image), image(image) {
    acc = 0;
}

void FractaleGLImage::idleFunc(){
    ++acc;

    //Do not change everytime to make something smoother
    if(acc == 1000){
       image->incN();
       updateView();
       acc = 0;
    }
}
