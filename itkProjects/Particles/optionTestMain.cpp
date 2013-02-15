//
//  optionTestMain.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/15/13.
//
//

#include <stdio.h>
#include <iostream>
#include "sstream"
#include "piOptions.h"

using namespace pi;
using namespace std;

int findLiteralType(string s) {
    stringstream ss(s);
    int v1;
    double v2;
    string v3;
    ss >> v1;
    if (ss.fail()) {
        ss >> v2;
        if (ss.fail()) {
            ss >> v3;
            if (v3 == "true" || v3 == "false") {
                return 4;
            }
            return 3;
        }
        return 2;
    }
    return 1;
}

int main(int argc, char* argv[]) {
    Options args;
    args.AppendString("args", "a");
    args.AppendString("hello", "b");
    args.AppendString("args", "c");
    args.Set("mother", 1);
    args.Set("father", 1.1);
    args.Set("brother", "Joseph");
    args.AppendInt("mine:", 1);
    args.AppendInt("mine:", 2);
    args.AppendInt("mine:", 3);
    args.AppendDouble("time:", 1.0);
    args.AppendDouble("time:", 2.3);
    args.AppendDouble("time:", 3.1);
    cout << args;

    stringstream os;
    os << args;

    Options args2;
    os >> args2;

    cout << "Save and Loaded: " << endl;
    cout << args2;
    cout << "Done" << endl;

    cout << findLiteralType("1.3f") << endl;
    cout << findLiteralType("1") << endl;
    cout << findLiteralType("hello") << endl;
    cout << findLiteralType("true") << endl;

    stringstream os2("hello int 1\n\n\nbyebye string joohwi\n");
    Options args3;
    os2 >> args3;
    cout << args3;
    cout << "end" << endl;
}