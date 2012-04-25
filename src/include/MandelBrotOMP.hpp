#include "Fractales.hpp"

#define THREADS 12

class MandelBrotImageOMP : public FractaleImage {
    public:
	MandelBrotImageOMP(int m, int n, DomaineMaths domain);

    protected:
	void refreshAll(const DomaineMaths& domainNew);

    private:
	float mandelbrot(float x, float y);
};
