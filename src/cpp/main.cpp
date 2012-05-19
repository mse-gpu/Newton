#include <iostream>
#include <limits.h>
#include <cmath>

#include "omp.h"
#include "ChronoOMPs.h"

#include "Fractales.hpp"
#include "NewtonSequential.hpp"
#include "NewtonOMP.hpp"

//The Newton Fractales Launchers
int launchNewton();
int launchNewtonOMP();

//The benchmark
int bench();

int main(void){
    //return launchNewton();
    //return launchNewtonOMP();
    return bench();
}

int launchNewton(){
    std::cout << "Launch Newton" << std::endl;

    char** argv = NULL;
    GLUTWindowManagers::init(0, argv);

    int w = 800;
    int h = 800;

    DomaineMaths domain(-w / 2, - h / 2, w, h);

    NewtonImageSequential* functionalImage = new NewtonImageSequential(w, h, domain);
    FractaleGLImage* functionSelections = new FractaleGLImage(functionalImage);

    GLUTWindowManagers* windowManager = GLUTWindowManagers::getInstance();
    windowManager->createWindow(functionSelections);
    windowManager->runALL(); //This call is blocking

    return 0;
}

int launchNewtonOMP(){
    omp_set_num_threads(THREADS);

    std::cout << "Launch Newton with OMP" << std::endl;

    char** argv = NULL;
    GLUTWindowManagers::init(0, argv);

    int w = 800;
    int h = 800;

    DomaineMaths domain(-w / 2, - h / 2, w, h);

    NewtonImageOMP* functionalImage = new NewtonImageOMP(w, h, domain);
    FractaleGLImage* functionSelections = new FractaleGLImage(functionalImage);

    GLUTWindowManagers* windowManager = GLUTWindowManagers::getInstance();
    windowManager->createWindow(functionSelections);
    windowManager->runALL(); //This call is blocking

    return 0;
}

#define DIM_H 10000
#define DIM_W 10000
#define TIMES 20

#define THREADS 24

namespace newton_bench {
#include "Newton.hpp"
}

struct rgba {
	int r;
	int g;
	int b;
	int a;
};

void setFloatRGBA(rgba* image, int i, int j, int r, int g, int b, int a){
    image[i * (DIM_H) + j].r = r;
    image[i * (DIM_H) + j].g = g;
    image[i * (DIM_H) + j].b = b;
    image[i * (DIM_H) + j].a = a;
}

void benchParallelNewton(rgba* image){
    omp_set_num_threads(THREADS);

    float xMin = -1.7;
    float xMax = +1.7;
    float yMin = -1.1;
    float yMax = +1.1;

    DomaineMaths domain(xMin, yMin, xMax - xMin, yMax - yMin);

    float dx = (float) (domain.dx / (float) DIM_W);
    float dy = (float) (domain.dy / (float) DIM_H);

    #pragma omp parallel
    {
    	int tid = omp_get_thread_num();
    	int i = tid + 1;

    	float y = domain.y0 + tid * dy;

    	while(i <= DIM_H){
    	    float x = domain.x0;

    	    for(int j = 1; j <= DIM_W; ++j){
		int color = newton_bench::real_newton(x, y);

		if(color == 0){
		    setFloatRGBA(image, i, j, 0, 0, 0, 0);
		} else if(color == 1){
		    setFloatRGBA(image, i, j, 1, 0, 0, 0);
		} else if(color == 2){
		    setFloatRGBA(image, i, j, 0, 1, 0, 0);
		} else if(color == 3){
		    setFloatRGBA(image, i, j, 0, 0, 1, 0);
		}

    		x += dx;
    	    }

    	    y += THREADS * dy;

    	    i += THREADS;
    	}
    }
}

void benchSequentialNewton(rgba* image){
    float xMin = -1.7;
    float xMax = +1.7;
    float yMin = -1.1;
    float yMax = +1.1;

    DomaineMaths domain(xMin, yMin, xMax - xMin, yMax - yMin);

    float dx = (float) (domain.dx / (float) DIM_W);
    float dy = (float) (domain.dy / (float) DIM_H);
    float y = domain.y0;

    float acc = 0;

    for(int i = 1; i <= DIM_H; ++i){
    	float x = domain.x0;

    	for(int j = 1; j <= DIM_W; ++j){
	    int color = newton_bench::real_newton(x, y);

	    if(color == 0){
		setFloatRGBA(image, i, j, 0, 0, 0, 0);
	    } else if(color == 1){
		setFloatRGBA(image, i, j, 1, 0, 0, 0);
	    } else if(color == 2){
		setFloatRGBA(image, i, j, 0, 1, 0, 0);
	    } else if(color == 3){
		setFloatRGBA(image, i, j, 0, 0, 1, 0);
	    }

    	    x += dx;
    	}

    	y += dy;
    }
}

void benchNewton(){
    std::cout << "Launch the Newton benchmark" << std::endl;

    rgba* image = (rgba*) calloc(sizeof(rgba),  (DIM_H + 1) * (DIM_W + 1));

    ChronoOMPs chronos;
    chronos.start();

    for(int i = 0; i < TIMES; ++i){
	benchSequentialNewton(image);
    }

    double timeSequential = chronos.timeElapse();
    std::cout << "Sequential Total (" << TIMES << " times) " << timeSequential << "s" << std::endl;
    std::cout << "Sequential Mean  (" << TIMES << " times) " << (timeSequential / TIMES) << "s" << std::endl;

    chronos.start();

    for(int i = 0; i < TIMES; ++i){
	benchParallelNewton(image);
    }

    double timeParallel = chronos.timeElapse();
    std::cout << "OMP Total (" << TIMES << " times) " << timeParallel << "s" << std::endl;
    std::cout << "OMP Mean  (" << TIMES << " times) " << (timeParallel / TIMES) << "s" << std::endl;

    std::cout << "Factor=" << (timeSequential / timeParallel) << std::endl;

    free(image);
}

int bench(){
    benchNewton();

    return 0;
}
