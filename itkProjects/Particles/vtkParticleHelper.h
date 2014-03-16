#ifndef VTKPARTICLEHELPER_H
#define VTKPARTICLEHELPER_H

#include "piParticleCore.h"
#include "piOptions.h"

#include "vtkPolyData.h"

class vtkParticleHelper
{
public:
    vtkParticleHelper();

    /// main test code
    bool main(pi::Options& opts, pi::StringVector& args);

    /// split the given polygon using loop subdivision filter
    void splitParticles(pi::ParticleSubject& subj);

private:
    void testSplitParticles();

};

#endif // VTKPARTICLEHELPER_H
