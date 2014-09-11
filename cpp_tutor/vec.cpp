#include <vector>
#include <iostream>

using namespace std;

int main1() {
	std::vector<int> arr(10);
	for (int j = 0; j < 100; j++) {
		arr[10+10*j] = 1;
	}
	for (int i = 0; i < 10+1000+1; i++) {
		cout << arr[10+1000] << ",";
	}
	cout << endl;
	return 0;
}

int main2() {
	std::vector<int> arr(3);
	for (int j = 0; j < arr.size(); j++) {
		arr[j] = j+1;
	}
	int* p = &arr[0];
	while (p != &arr[3]) {
		cout << *p << ",";
		p++;
	}
	cout << endl;
	std::vector<int>::iterator iter = arr.begin();
	while (iter != arr.end()) {
		cout << *iter << ",";
		iter++;
	}
	cout << endl;
	return 0;
}

#include <algorithm>
void print(int x) {
	cout << x << ",";
}
int main3() {
	vector<int> a(3);
	std::fill(a.begin(), a.end(), 100);
	std::for_each(a.begin(), a.end(), print);
	cout << endl;
	return 0;
}


struct Odd2Counter {
	int num;
	Odd2Counter(int s):num(s) {
	}

	int operator()() {
		num += 1;
		return num*num;
	}
};

int main4() {
	vector<int> a(3);
	Odd2Counter counter(0);
	std::generate(a.begin(), a.end(), counter);
	std::for_each(a.begin(), a.end(), print);
	cout << endl;
	return 0;
}


bool is_odd(int i) {
	return i%2 == 0;
}
int main5() {
	vector<int> a(4);
	Odd2Counter counter(0);
	std::generate(a.begin(), a.end(), counter);
	if (std::any_of(a.begin(), a.end(), is_odd)) {
		cout << "an odd number detected" << endl;
	}

	if (std::all_of(a.begin(), a.end(), is_odd)) {
		cout << "all odd numbers" << endl;
	}

	cout << "# odds: " << std::count_if(a.begin(), a.end(), is_odd) << endl;
	return 0;
}


int add(int x, int y) {
	return x + y;
}
int sqr(int x) {
	return x*x;
}
int main() {
	vector<int> a(5), b(5);
	std::fill(a.begin(), a.end(), 1);
	std::transform(a.begin()+1, a.end(), a.begin(), a.begin()+1, add);
	std::for_each(a.begin(), a.end(), print);
	cout << endl;
	std::transform(a.begin(), a.end(), a.begin(), sqr);
	std::for_each(a.begin(), a.end(), print);
	cout << endl;
	return 0;
}
