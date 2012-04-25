#include <iostream>
#include <cmath>

#include "omp.h"

#include "JuliaSequential.hpp"

JuliaImageSequential::JuliaImageSequential(int m, int n, DomaineMaths domain, float cReal, float cImag) : FractaleImage(m,n,domain,10), cReal(cReal), cImag(cImag) {
    //Nothing to init
}

void JuliaImageSequential::refreshAll(const DomaineMaths& domainNew){
    int w = getW();
    int h = getH();

    float dx = (float) (domainNew.dx / (float) w);
    float dy = (float) (domainNew.dy / (float) h);
    float y = domainNew.y0;

    for(int i = 1; i <= h; ++i){
	float x = domainNew.x0;

	for(int j = 1; j <= w; ++j){
	    float h = julia(x, y);

	    //setFloatRGBA(i, j, h, h, h);
	    if(h == 0){
		setHSB(i, j, 0, 0, 0);
	    } else {
		setHSB(i, j, h, 1.0, 1.0);
	    }

	    x += dx;
	}

	y += dy;
    }
}

float JuliaImageSequential::julia(float x, float y){
    float real = x;
    float imag = y;

    float n = 0;
    float norm;

    do{
	float tmpReal = real;
	real = real * real - imag * imag + cReal;
	imag = tmpReal * imag + imag * tmpReal + cImag;

	++n;

	norm = sqrt(real * real + imag * imag);
    } while (norm <= 2.0 && n < N);

    return n == N ? 0 : (n / (float) N);
}
