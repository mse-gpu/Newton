#ifndef FRACTALES
#define FRACTALES

#include "DomaineMaths.h"
#include "GLUTWindowManagers.h"
#include "ImageFonctionelSelectionMOOs.h"
#include "GLImageFonctionelSelections.h"

class FractaleImage : public ImageFonctionelSelectionMOOs {
    public:
	FractaleImage(int m, int n, DomaineMaths domain, int N);

	void incN();

    protected:
	void onDomaineChangePerformed(const DomaineMaths& domainNew);
	virtual void refreshAll(const DomaineMaths& domainNew) = 0;

	int N;
};

class FractaleGLImage : public GLImageFonctionelSelections {
    public:
	FractaleGLImage(FractaleImage* image);

    private:
	float acc;

	void idleFunc();

	FractaleImage* image;
};

#endif
