#include <iostream>

using namespace std;


int main() {
	char** buf = new (std::nothrow) char*[1];
	buf[0] = (char*) 0x1;
	buf = (char**) 0x1;
	/*
	buf = new char[10];

	for (int j = 0; j < 10; j++) {
		cout << buf[j] << ",";
	}
	cout << endl;
	*/
	delete[] buf;
}
