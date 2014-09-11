#include <exception>

class Hello {
	public:
		Hello() {
			throw new std::exception();
		}
};


int main() {
	try {
		Hello hello;
	} catch (std::exception* e) {
	}
}
