#include <iostream>
#include <cmath>

#include "omp.h"

#include "NewtonOMP.hpp"

//TODO That's ugly...
namespace newton_omp {
#include "Newton.hpp"
}

NewtonImageOMP::NewtonImageOMP(int m, int n, DomaineMaths domain) : FractaleImage(m,n,domain) {
    refreshAll(domain);
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
		int color = newton_omp::real_newton(x, y);

		if(color == 0){
		    setFloatRGBA(i, j, 0, 0, 0);
		} else if(color == 1){
		    setFloatRGBA(i, j, 1, 0, 0);
		} else if(color == 2){
		    setFloatRGBA(i, j, 0, 1, 0);
		} else if(color == 3){
		    setFloatRGBA(i, j, 0, 0, 1);
		}

		x += dx;
	    }

	    y += THREADS * dy;

	    i += THREADS;
	}
    }
}

