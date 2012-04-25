#include <iostream>
#include <limits.h>
#include <cmath>

#include "omp.h"
#include "ChronoOMPs.h"

#include "Fractales.hpp"
#include "JuliaSequential.hpp"
#include "JuliaOMP.hpp"
#include "MandelBrotSequential.hpp"
#include "MandelBrotOMP.hpp"

//The MandelBrot Fractales Launchers
int launchMandelbrot();
int launchMandelbrotOMP();

//The Julia Fractales Launchers
int launchJulia();
int launchJuliaOMP();

//The benchmark
int bench();

int main(void){
    //return launchMandelbrot();
    //return launchMandelbrotOMP();
    //return launchJulia();
    //return launchJuliaOMP();
    return bench();
}

int launchMandelbrot(){
    std::cout << "Launch the application" << std::endl;

    char** argv = NULL;
    GLUTWindowManagers::init(0, argv);

    float xMin = -1.3968;
    float xMax = -1.3578;
    float yMin = -0.03362;
    float yMax = 0.0013973;

    DomaineMaths domain(xMin, yMin, xMax - xMin, yMax - yMin);

    int w = 800;
    int h = 600;

    MandelBrotImageSequential* functionalImage = new MandelBrotImageSequential(w, h, domain);
    FractaleGLImage* functionSelections = new FractaleGLImage(functionalImage);

    GLUTWindowManagers* windowManager = GLUTWindowManagers::getInstance();
    windowManager->createWindow(functionSelections);
    windowManager->runALL(); //This call is blocking

    return 0;
}

int launchMandelbrotOMP(){
    omp_set_num_threads(THREADS);

    std::cout << "Launch the application with OMP" << std::endl;

    char** argv = NULL;
    GLUTWindowManagers::init(0, argv);

    float xMin = -1.3968;
    float xMax = -1.3578;
    float yMin = -0.03362;
    float yMax = 0.0013973;

    DomaineMaths domain(xMin, yMin, xMax - xMin, yMax - yMin);

    int w = 800;
    int h = 600;

    MandelBrotImageOMP* functionalImage = new MandelBrotImageOMP(w, h, domain);
    FractaleGLImage* functionSelections = new FractaleGLImage(functionalImage);

    GLUTWindowManagers* windowManager = GLUTWindowManagers::getInstance();
    windowManager->createWindow(functionSelections);
    windowManager->runALL(); //This call is blocking

    return 0;
}

int launchJulia(){
    std::cout << "Launch Julia" << std::endl;

    char** argv = NULL;
    GLUTWindowManagers::init(0, argv);

    float xMin = -1.7;
    float xMax = +1.7;
    float yMin = -1.1;
    float yMax = +1.1;

    DomaineMaths domain(xMin, yMin, xMax - xMin, yMax - yMin);

    int w = 800;
    int h = 600;

    JuliaImageSequential* functionalImage = new JuliaImageSequential(w, h, domain, -0.745, +0.1);
    FractaleGLImage* functionSelections = new FractaleGLImage(functionalImage);

    GLUTWindowManagers* windowManager = GLUTWindowManagers::getInstance();
    windowManager->createWindow(functionSelections);
    windowManager->runALL(); //This call is blocking

    return 0;
}

int launchJuliaOMP(){
    omp_set_num_threads(THREADS);

    std::cout << "Launch Julia with OMP" << std::endl;

    char** argv = NULL;
    GLUTWindowManagers::init(0, argv);

    float xMin = -1.7;
    float xMax = +1.7;
    float yMin = -1.1;
    float yMax = +1.1;

    DomaineMaths domain(xMin, yMin, xMax - xMin, yMax - yMin);

    int w = 800;
    int h = 600;

    JuliaImageOMP* functionalImage = new JuliaImageOMP(w, h, domain, -0.745, +0.1);
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

float mandelbrot(float x, float y){
    float imag = 0.0;
    float real = 0.0;

    float n = 0;
    float norm;

    do{
	float tmpReal = real;
	real = real * real - imag * imag + x;
	imag = tmpReal * imag + imag * tmpReal + y;

	++n;

	norm = sqrt(real * real + imag * imag);
    } while (norm <= 2.0 && n < N);

    return n == N ? 0 : (n / (float) N);
}

void benchSequentialMandelBrot(){
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
	    float h = mandelbrot(x, y);

	    acc += h;

	    x += dx;
	}

	y += dy;
    }
}

void benchParallelMandelBrot(){
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

        float acc = 0;

    	float y = domain.y0 + tid * dy;

    	while(i <= DIM_H){
    	    float x = domain.x0;

    	    for(int j = 1; j <= DIM_W; ++j){
    		float h = mandelbrot(x, y);

    		acc += h;

    		x += dx;
    	    }

    	    y += THREADS * dy;

    	    i += THREADS;
    	}
    }
}

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

void benchParallelJulia(){
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

void benchSequentialJulia(){
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

void benchJulia(){
    std::cout << "Launch the Julia benchmark" << std::endl;

    ChronoOMPs chronos;
    chronos.start();

    benchSequentialJulia();

    double timeSequential = chronos.timeElapse();
    std::cout << "Sequential version took " << timeSequential << "s" << std::endl;

    chronos.start();

    benchParallelJulia();

    double timeParallel = chronos.timeElapse();
    std::cout << "OMP version took " << timeParallel << "s" << std::endl;

    std::cout << "Factor=" << (timeSequential / timeParallel) << std::endl;
}

void benchMandelbrot(){
    std::cout << "Launch the MandelBrot benchmark" << std::endl;

    ChronoOMPs chronos;
    chronos.start();

    benchSequentialMandelBrot();

    double timeSequential = chronos.timeElapse();
    std::cout << "Sequential version took " << timeSequential << "s" << std::endl;

    chronos.start();

    benchParallelMandelBrot();

    double timeParallel = chronos.timeElapse();
    std::cout << "OMP version took " << timeParallel << "s" << std::endl;

    std::cout << "Factor=" << (timeSequential / timeParallel) << std::endl;
}

int bench(){
    benchMandelbrot();
    benchJulia();

    return 0;
}
