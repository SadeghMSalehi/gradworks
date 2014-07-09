#include "test1.h"
#include <iostream>


using namespace std;

class A {
public:
    int id;
    string name;
    double height;


private:
    friend class B;
    virtual void hello() {
        cout << "Hello" << endl;
    }

public:
    void hellox() {
        cout << "Hello A" << endl;
    }
    virtual void helloy() {
        cout << "Hello A~" << endl;
    }
};

class AA: public A {
private:
    virtual void hello() {
        cout << "Hello AA" << endl;
    }

public:
    void hellox() {
        cout << "hello AA" << endl;
    }
    virtual void helloy() {
        cout << "Hello AA~" << endl;
    }
};


class B {
public:
    void hello() {
//        AA a;
        A a;
        a.hello();
    }
};

void run_friend_test() {
    B b;
    b.hello();
}

void run_virtual_test() {
    A* a = new A();
    AA* aa = new AA();

    A* xa = (A*) aa;

    // call non-virtual hellox() of A
    a->hellox();
    // call non-virtual hellox() of A even though it is aa
    xa->hellox();
    // call non-virtual hellox() of AA
    aa->hellox();


    // call virtual helloy() of A
    a->helloy();
    // call virtual helloy() of AA (its original form)
    xa->helloy();
}


namespace animal {
    class Dog {
    protected:
        void say() {
            cout << "bow" << endl;
        }
    private:
        void name() {
            cout << "unknown" << endl;
        }
    };
    class BullDog: public Dog {
    public:
        void wow() {
            say();
//            name();
        }
    };
    class DogOwner {
    public:
        void say() {
            Dog d;
//            d.say();
            cout << "hello";
            
        }
    };
}

namespace botanics {

}

void run_scope_test() {
    animal::BullDog dog;
    dog.wow();
}

void run_init_test() {
//    A a = { 172, "name", 3.0 };
    
}

int main(int argc, char* argv[]) {
    run_friend_test();
    run_virtual_test();
    run_scope_test();
    run_init_test();
    return 0;
}