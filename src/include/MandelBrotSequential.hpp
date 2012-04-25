#include "Fractales.hpp"

class MandelBrotImageSequential : public FractaleImage {
    public:
	MandelBrotImageSequential(int m, int n, DomaineMaths domain);

    protected:
	void refreshAll(const DomaineMaths& domainNew);

    private:
	float mandelbrot(float x, float y);
};
