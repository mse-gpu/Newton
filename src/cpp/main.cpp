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
    return launchNewtonOMP();
    //return bench();
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

#define DIM_H 12000
#define DIM_W 16000

#define THREADS 12

#define N 52

static float cReal = -0.745;
static float cImag = +0.1;

float julia(float x, float y){
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

void benchParallelNewton(){
    omp_set_num_threads(THREADS);

    float xMin = -1.7;
    float xMax = +1.7;
    float yMin = -1.1;
    float yMax = +1.1;

    DomaineMaths domain(xMin, yMin, xMax - xMin, yMax - yMin);

    float dx = (float) (domain.dx / (float) DIM_W);
    float dy = (float) (domain.dy / (float) DIM_H);
    float acc = 0;

    #pragma omp parallel
    {
    	int tid = omp_get_thread_num();
    	int i = tid + 1;

    	float y = domain.y0 + tid * dy;

    	while(i <= DIM_H){
    	    float x = domain.x0;

    	    for(int j = 1; j <= DIM_W; ++j){
    		float h = julia(x, y);

    		acc += h;

    		x += dx;
    	    }

    	    y += THREADS * dy;

    	    i += THREADS;
    	}
    }
}

void benchSequentialNewton(){
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
    	    float h = julia(x, y);

    	    acc += h;

    	    x += dx;
    	}

    	y += dy;
    }
}

void benchNewton(){
    std::cout << "Launch the Newton benchmark" << std::endl;

    ChronoOMPs chronos;
    chronos.start();

    benchSequentialNewton();

    double timeSequential = chronos.timeElapse();
    std::cout << "Sequential version took " << timeSequential << "s" << std::endl;

    chronos.start();

    benchParallelNewton();

    double timeParallel = chronos.timeElapse();
    std::cout << "OMP version took " << timeParallel << "s" << std::endl;

    std::cout << "Factor=" << (timeSequential / timeParallel) << std::endl;
}

int bench(){
    benchNewton();

    return 0;
}
