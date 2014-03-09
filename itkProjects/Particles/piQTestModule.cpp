//
//  piQTestModule.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/4/13.
//
//

#include "piqTestModule.h"
#include "piParticleTrace.h"
#include "piFitCurve.h"
#include "piBSplineBasis.h"
#include "piContourSystem.h"
#include "piImageIO.h"
#include <itkDomainThreader.h>
#include <itkThreadedImageRegionPartitioner.h>
#include <itkImageRegionConstIterator.h>

using namespace std;
using namespace pi;
using namespace itk;

namespace piq {
    class DomainCounter {
    public:
        
    };

    class DomainCounterThreader: public DomainThreader<ThreadedImageRegionPartitioner<3>, DomainCounter> {
    public:
        typedef DomainCounterThreader Self;
        typedef DomainThreader Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        typedef Superclass::DomainType DomainType;

        itkNewMacro(DomainCounterThreader);
        itkTypeMacro(DomainCounterThreader, DomainThreader);

    protected:
        DomainCounterThreader() {}
        virtual ~DomainCounterThreader() {}
        void BeforeThreadedExecution() {
            cout << "Before Execution" << endl;
        }
        void ThreadedExecution(const DomainType& subDomain, const ThreadIdType threadId) {
            cout << threadId << endl;
        }
        void AfterThreadedExecution() {
            cout << "After Execution" << endl;
        }

    private:
        DomainCounterThreader(const DomainCounterThreader&);
        void operator=(const DomainCounterThreader&);
    };

    TestModule::TestModule(QObject* parent): QObject(parent) {

    }

    TestModule::~TestModule() {
        
    }

    void TestModule::ThreadTest() {
        typedef ThreadedImageRegionPartitioner<3> PartitionType;
        PartitionType::Pointer part = PartitionType::New();
        PartitionType::DomainType region;
        for (int i = 0; i < 3; i++) {
            region.SetIndex(i, 0);
            region.SetSize(i, 100);
        }

        PartitionType::DomainType subRegion;
        for (int i = 0; i < 10; i++) {
            part->PartitionDomain(i, 10, region, subRegion);
//            cout << (i+1) << "/10:" << subRegion << endl;
        }

        DomainCounterThreader::Pointer domainCounter = DomainCounterThreader::New();
        domainCounter->SetMaximumNumberOfThreads(2);
        domainCounter->Execute(NULL, region);
    }

    void TestModule::ContourTest() {
        ifstream is("/NIRAL/work/joohwi/nadia/C31_E04_slice128_curve.txt");
        ParticleVector particles;
        ParticleTrace::Read(is, particles);

        CurveFitting fit;
        fit.FitCurve(particles);
        ParticleVector& result = fit.GetResult();

        ImageIO<RealSlice> io;
        RealSlice::Pointer image1 = io.ReadCastedImage("/NIRAL/work/joohwi/nadia/SliceImages/C31_E04_slice_128.nii.gz");
        RealSlice::Pointer image2 = io.ReadCastedImage("/NIRAL/work/joohwi/nadia/SliceImages/C31_E04_slice_129.nii.gz");

        RealSliceVector images;
        images.push_back(image1);
        images.push_back(image2);

        ContourSystem system;
        system.SetSlices(images);
        system.SetInitialParticles(result);
        system.SetAttributeDimensions(10, 20);
        system.AllocateAttributeBuffer();
        system.Track(1);
    }

    void TestModule::BSplineTest() {
        ifstream is("/NIRAL/work/joohwi/nadia/C31_E04_slice128_curve.txt");
        ParticleVector particles, controls, fittingResult;
        ParticleTrace::Read(is, particles);
        cout << "Read " << particles.size() << " particles ..." << endl;

        const int nControls = 20;
        BSplineBasis::CubicContourFitting(particles, nControls, controls);
        for (int i = 0; i < nControls; i++) {
            cout << controls[i] << endl;
        }

        std::vector<float> params(50);
        for (int i = 0; i < 50; i++) {
            params[i] = i / 50.0;
        }

        BSplineBasis::CubicContour(controls, params, fittingResult);

        cout << "Contour:" << endl;
        for (int i = 0; i < fittingResult.size(); i++) {
            cout << fittingResult[i] << endl;
        }
    }

    void TestModule::BSplineBasisTest() {
        const int k = 3;
        const int np = 3;
        const int nk = np + k + 1;
        BSplineBasis::TReal knots[nk] = { 0, 0, 0, 1, 2, 2, 2 };
        pi::BSplineBasis b1(knots, nk, 0, k);
        pi::BSplineBasis b2(knots, nk, 1, k);
        pi::BSplineBasis b3(knots, nk, 2, k);
        pi::BSplineBasis b4(knots, nk, 3, k);
        for (float t = 0; t < 2.1; t += 0.1) {
            cout << t << "\t" <<  b1(t) << "\t" << b2(t) << "\t" << b3(t) << "\t" << b4(t) << "\t" << endl;
        }
    }
    
    void TestModule::FitTest(pi::StringVector &args) {
        ifstream is(args[0].c_str());
        ParticleVector particles;
        ParticleTrace::Read(is, particles);
        cout << "Read " << particles.size() << " particles ..." << endl;

        CurveFitting fit;
        fit.FitCurve(particles);
        ParticleVector& result = fit.GetResult();

        ofstream os("/tmpfs/curvefit.txt");
        ParticleTrace::Write(os, result);
        os.close();
    }

    void TestModule::Run(Options *opts, StringVector &args) {
        if (opts->GetBool("--fitTest")) {
            FitTest(args);
        } else if (opts->GetBool("--bsplineBasisTest")) {
            BSplineTest();
        } else if (opts->GetBool("--contourTest")) {
            ContourTest();
        } else if (opts->GetBool("--threadTest")) {
            ThreadTest();
        }
    }
}