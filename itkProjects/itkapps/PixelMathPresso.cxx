#include "pthread.h"

#include "itkImageCommon.h"
typedef itk::Image<unsigned int,3> ImageType;

#include "itkUnaryFunctorImageFilter.h"

#include "boost/program_options.hpp"
namespace po = boost::program_options;

#include "string"
#include "vector"
#include "iostream"
using namespace std;

void print_help_message(char* argv[]) {
	printf("usage: %s [options] files \n\
	--eq,-e 	equation for PixelMath \n\
	--output,-o	output image file \n\
	\n\
	In the equation, pixels are represented as variable x,y,z (only x for now).\n\
	Take the variable as pixel input and the equation result will be the pixel output.\n\
	\n\
	Examples)\n\
		+1 for each pixel >> PixelMath -e \"x+1\" -o plus1.nrrd [yourinput.nrrd]\n\
		label swap 1 to 2; 2 to 1 >> PixelMath -e \"(x==1) ? 2 : ((x==2) ? 1 : x)\"\n\
", argv[0]);
}

#include <MathPresso/MathPresso.h>

static double plus1(double x) {
	return x+1;
}

void run(vector<string> &inputfiles, string &outfile, string &eq) {
	struct Variables {
		MathPresso::mreal_t x;
		MathPresso::mreal_t y;
		MathPresso::mreal_t z;
	};

  MathPresso::Context ctx;
  MathPresso::Expression e;

  Variables variables;
  variables.x = 0.0;
  variables.y = 0.0;
  variables.z = 0.0;

  ctx.addEnvironment(MathPresso::MENVIRONMENT_ALL);
  ctx.addVariable("x", MATHPRESSO_OFFSET(Variables, x));
  ctx.addVariable("y", MATHPRESSO_OFFSET(Variables, y));
  ctx.addVariable("z", MATHPRESSO_OFFSET(Variables, z));

	const int N9 = 1e9;
	const int N6 = 1e6;
	double *data = new double[N9];
	if (data == NULL) {
		cout << "out of memory" << endl;
		exit(0);
	}
	MathPresso::mresult_t result = e.create(ctx, eq.c_str());
	if (result != MathPresso::MRESULT_OK || result == MathPresso::MRESULT_NO_EXPRESSION) {
		fprintf(stderr, "Error compiling expression:\n%s\n", eq.c_str());
		return;
	} 

	for (int i = 0; i <N9; i++) {
		data[i] = i;
		variables.x = data[i];
//		data[i] = data[i] + 1;
		data[i] = e.evaluate(&variables);
//		data[i] = plus1(variables.x);
	}
}

int main(int argc, char* argv[]) {
	string outfile, eq;
	po::options_description desc("Options");
	desc.add_options()
		("help,h", "produce help messages")
		("eq,e", po::value<string>(&eq), "equation for pixels ex) x+1")		
		("output,o", po::value<string>(&outfile)->default_value("result.nrrd"), "output filename")
	;

	po::options_description args("Arguments");
	args.add_options()
		("files", po::value< vector<string> >(), "files... ")
	;

	po::positional_options_description p;
	p.add("files", -1);

	po::options_description all_opts("All options");
	all_opts.add(args).add(desc);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all_opts).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("files") < 1 || eq == "" || vm.count("help") > 0) {
		print_help_message(argv);
		return 0;
	}

	vector<string> inputFiles = vm["files"].as< vector<string> >();
	run(inputFiles, outfile, eq);
}
