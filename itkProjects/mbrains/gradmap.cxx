#include "piImageDef.h"
#include "piImageProcessing.h"
#include "iostream"
#include "itkImageIO.h"

using namespace std;
using namespace pi;
using namespace itkcmds;

int main(int argc, char* argv[]) {
    Options args;
    StringVector args = args.ParseArgs(argc, argv);
    ImageProcessing proc;

    itkImageIO<RealImage> real;

}