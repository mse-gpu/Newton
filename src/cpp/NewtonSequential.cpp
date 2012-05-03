#include <iostream>
#include <cmath>

#include "omp.h"

#include "NewtonSequential.hpp"

//TODO That's ugly...
namespace newton_seq {
#include "Newton.hpp"
}

NewtonImageSequential::NewtonImageSequential(int m, int n, DomaineMaths domain) : FractaleImage(m,n,domain){
    refreshAll(domain);
}

void NewtonImageSequential::refreshAll(const DomaineMaths& domainNew){
    int w = getW();
    int h = getH();

    float dx = (float) (domainNew.dx / (float) w);
    float dy = (float) (domainNew.dy / (float) h);
    float y = domainNew.y0;

    for(int i = 1; i <= h; ++i){
	float x = domainNew.x0;

	for(int j = 1; j <= w; ++j){
	    int color = newton_seq::real_newton(x, y);

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

	y += dy;
    }
}
