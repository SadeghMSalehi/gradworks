#include "kfem.h"
#include "laplaceEx.h"

using namespace std;
using namespace dealii;

int main(int argc, char* argv[]) {
	LaplaceEx lex;
	lex.run();
	return 0;
}
