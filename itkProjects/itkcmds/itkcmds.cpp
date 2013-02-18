#include "iostream"


#include "piOptions.h"
#include "itkImageIO.h"

using namespace pi;
int itkNormalizedCorrelationImageMetricTest(int, char* []);

int main(int argc, char* argv[]) {
    CSimpleOpt::SOption spec[] = {
        { 10, "--ncctest", SO_NONE },
        SO_END_OF_OPTIONS
    };
    
    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, spec);

    if (parser.GetBool("--ncctest")) {
        itkNormalizedCorrelationImageMetricTest(0,NULL);
        return 0;
    }
    
}