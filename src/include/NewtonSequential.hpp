#include "Fractales.hpp"

class NewtonImageSequential : public FractaleImage {
    public:
	NewtonImageSequential(int m, int n, DomaineMaths domain);

    protected:
	void refreshAll(const DomaineMaths& domainNew);
};
