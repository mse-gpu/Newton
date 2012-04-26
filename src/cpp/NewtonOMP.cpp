#include <iostream>
#include <cmath>

#include "omp.h"

#include "NewtonOMP.hpp"

NewtonImageOMP::NewtonImageOMP(int m, int n, DomaineMaths domain, float cReal, float cImag) : FractaleImage(m,n,domain), cReal(cReal), cImag(cImag) {
    //Nothing to init
}

void NewtonImageOMP::refreshAll(const DomaineMaths& domainNew){
    int w = getW();
    int h = getH();

    float dx = (float) (domainNew.dx / (float) w);
    float dy = (float) (domainNew.dy / (float) h);

    #pragma omp parallel
    {
	int tid = omp_get_thread_num();
	int i = tid + 1;

	float y = domainNew.y0 + tid * dy;

	while(i <= h){
	    float x = domainNew.x0;

	    for(int j = 1; j <= w; ++j){
		float h = newton(x, y);

		//setFloatRGBA(i, j, h, h, h);
		if(h == 0){
		    setHSB(i, j, 0, 0, 0);
		} else {
		    setHSB(i, j, h, 1.0, 1.0);
		}

		x += dx;
	    }

	    y += THREADS * dy;

	    i += THREADS;
	}
    }
}

float NewtonImageOMP::newton(float x, float y){
    return 0;
}
