#include "itkImageIO.h"
#include "itkImage.h"
#include "iostream"

typedef itk::Image<double,3> Image;

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cout << "usage: " << argv[0] << " input.lst image.nrrd output.txt" << endl;
        return 0;
    }

    itkcmds::itkImageIO<Image> io;
    Image::Pointer img = io.ReadImageT(argv[2]);

    ifstream in(argv[1]);
    if (!in.is_open()) {
        return 0;
    }

    ofstream out(argv[3]);
    while (in.good()) {
        double x, y, z;
        in >> x >> y >> z;

        // assume physical point
        Image::PointType inPoint;
        inPoint[0] = x;
        inPoint[1] = y;
        inPoint[2] = z;

        Image::IndexType outIdx;

        img->TransformPhysicalPointToIndex(inPoint, outIdx);
        double pixel = img->GetPixel(outIdx);
        out << pixel;

        inPoint[0] = -x;
        inPoint[1] = -y;
        inPoint[2] = z;

        img->TransformPhysicalPointToIndex(inPoint, outIdx);
        pixel = img->GetPixel(outIdx);
        out << " " << pixel;

        out << endl;
    }

    out.close();
    in.close();
}
