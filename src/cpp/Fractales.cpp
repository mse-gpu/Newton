#include <iostream>
#include <cmath>

#include "Fractales.hpp"

FractaleImage::FractaleImage(int m, int n, DomaineMaths domain) : ImageFonctionelSelectionMOOs(m,n,domain){
    //Nothing to init
}

void FractaleImage::onDomaineChangePerformed(const DomaineMaths& domainNew){
    //Repaint everything
    refreshAll(domainNew);
}

FractaleGLImage::FractaleGLImage(FractaleImage* image) : GLImageFonctionelSelections(image), image(image) {
    //Nothing to init
}

void FractaleGLImage::idleFunc(){
    //Nothing
    //TODO Remove the function
}
