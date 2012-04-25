#include "Fractales.hpp"

#define THREADS 12

class NewtonImageOMP : public FractaleImage {
    public:
	NewtonImageOMP(int m, int n, DomaineMaths domain, float cReal, float cImag);

    protected:
	void refreshAll(const DomaineMaths& domainNew);

    private:
	float newton(float x, float y);

	const float cReal;
	const float cImag;
};
