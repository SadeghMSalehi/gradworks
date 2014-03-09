#include <iostream>
#include "itkImageIOBase.h"
#include "itkImage.h"

#include "itkImageIO.h"

typedef itk::Image<int, 3> IntImage;

using namespace std;
using namespace itkcmds;
using namespace itk;


#define forD(v) for (int v = 0; v < numDimensions; v++)
static void printImageInformation(ImageIOBase::Pointer imageBase) {
    cout << imageBase->GetNumberOfComponents() << "\t";
    int numDimensions = imageBase->GetNumberOfDimensions();
    cout << numDimensions << "\t";
    forD(i) {
        cout << imageBase->GetOrigin(i) << " ";
    }
    cout << "\t";
    forD(i) {
        cout << imageBase->GetSpacing(i) << " ";
    }
    cout << "\t";
    cout << imageBase->GetPixelTypeAsString(imageBase->GetPixelType());
    cout << "\t";
    cout << imageBase->GetComponentTypeAsString(imageBase->GetComponentType());
    // cout << "[" << imageBase->GetComponentType() << "]";
    cout << "\t";
    forD(i) {
        if (i > 0) {
            cout << ";";
        }
        forD(j) {
            if (j > 0) {
                cout << " ";
            }
            cout << imageBase->GetDirection(i)[j];
        }
    }
    cout << "\t";
    forD(i) {
        if (i > 0) {
            cout << " ";
        }
        cout << imageBase->GetDimensions(i);
    }
    cout << "\t";
    forD(i) {
        if (i > 0) {
            cout << " ";
        }
        cout << (imageBase->GetDimensions(i)*imageBase->GetSpacing(i));
    }
    cout << "\n";
    return;
}

int mainImageInfo(int argc, const char* argv[]) {
    cout << "FileName\t#components\t#dimensions\torigin\tspacing\tpixel\tcomponent\tdirection\tsize\tFOV" << endl;
    for (int i = 1; i < argc; i++) {
        itkImageIO<IntImage> io;
        ImageIOBase::Pointer imageBase = io.ReadImageInfo(argv[i]);
        cout << argv[i] << "\t";
        printImageInformation(imageBase);
    }
	return 0;
}
