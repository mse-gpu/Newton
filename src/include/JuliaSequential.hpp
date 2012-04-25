#include "Fractales.hpp"

class JuliaImageSequential : public FractaleImage {
    public:
	JuliaImageSequential(int m, int n, DomaineMaths domain, float cReal, float cImag);

    protected:
	void refreshAll(const DomaineMaths& domainNew);

    private:
	float julia(float x, float y);

	const float cReal;
	const float cImag;
};
