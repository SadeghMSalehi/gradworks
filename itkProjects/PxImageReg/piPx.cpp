//
//  piPx.cpp
//  PxImageReg
//
//  Created by Joohwi Lee on 11/5/13.
//
//

#include "piPx.h"

namespace pi {
    using namespace std;
    
    // utility operator overloading
    ostream& operator<<(ostream& os, const Px& par) {
        fordim(k) { os << par.x[k] << " "; }
        return os;
    }
    ostream& operator<<(ostream& os, const Px::Vector& par) {
        Px::Vector::const_iterator p = par.begin();

        if (p == par.end()) {
            return os;
        }

        os << *p;
        for (p++; p != par.end(); p++) {
            os << endl << *p;
        }

        return os;
    }

}