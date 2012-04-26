#include "Fractales.hpp"

#define THREADS 12

class NewtonImageOMP : public FractaleImage {
    public:
	NewtonImageOMP(int m, int n, DomaineMaths domain);

    protected:
	void refreshAll(const DomaineMaths& domainNew);

    private:
	int newton(float x, float y);
};
