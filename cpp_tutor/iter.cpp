#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int value() {
	static int v = 0;
	return v++;
}

void print(int a) {
	cout << a << ",";
}

int main() {
	vector<int> a(3);
	generate(a.begin(), a.end(), value);
	vector<int>::reverse_iterator iter = a.rbegin();
	for (; iter != a.rend(); iter++) {
		cout << *iter;
	}
	cout << endl;
	for_each(a.begin(), a.end(), print);
	cout << endl;
	return 0;
}

