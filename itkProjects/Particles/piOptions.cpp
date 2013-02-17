//
//  myOptions.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "piOptions.h"
#include "sstream"
#include "iostream"
#include "fstream"

using namespace std;

namespace pi {
    void Options::SetBool(std::string name, bool value) {
        std::pair<BoolMap::iterator, bool> result = _boolMap.insert(BoolPair(name, value));
        if (!result.second) {
            result.first->second = value;
        }
    }
    
    void Options::SetInt(std::string name, int value) {
        std::pair<IntMap::iterator, bool> result = _intMap.insert(IntPair(name, value));
        if (!result.second) {
            result.first->second = value;
        }
    }

    void Options::SetDouble(std::string name, double value) {
        std::pair<DoubleMap::iterator, bool> result = _doubleMap.insert(DoublePair(name, value));
        if (!result.second) {
            result.first->second = value;
        }    }

    void Options::SetString(std::string name, std::string value) {
        std::pair<StringMap::iterator, bool> result = _stringMap.insert(StringPair(name, value));
        if (!result.second) {
            result.first->second = value;
        }
    }

    void Options::AppendString(std::string name, std::string value) {
        _stringVectorMap[name].push_back(value);
    }

    void Options::AppendDouble(std::string name, double value) {
        _doubleVectorMap[name].push_back(value);
    }

    void Options::AppendInt(std::string name, int value) {
        _intVectorMap[name].push_back(value);
    }

    bool Options::GetBool(std::string name, bool def) {
        if (_boolMap.find(name) == _boolMap.end()) {
            return def;
        } else {
            return _boolMap.at(name);
        }
    }


    int Options::GetInt(std::string name, int def) {
        if (_intMap.find(name) == _intMap.end()) {
            return def;
        } else {
            return _intMap.at(name);
        }
    }

    int Options::GetStringAsInt(std::string name, int def) {
        return def;
    }

    double Options::GetDouble(std::string name, double def) {
        if (_doubleMap.find(name) == _doubleMap.end()) {
            return def;
        } else {
            return _doubleMap.at(name);
        }
    }

    double Options::GetStringAsDouble(std::string name, double def) {
        return def;
    }

    std::string Options::GetString(std::string name, std::string def) {
        if (_stringMap.find(name) == _stringMap.end()) {
            return def;
        } else {
            return _stringMap.at(name);
        }
    }

    StringVector& Options::GetStringVector(std::string name) {
        return _stringVectorMap[name];
    }

    std::string Options::GetStringVectorValue(std::string name, int n, std::string def) {
        StringVector& vector = _stringVectorMap[name];
        if (vector.size() > n) {
            return vector[n];
        }
        return def;
    }

    DoubleVector& Options::GetDoubleVector(std::string name) {
        return _doubleVectorMap[name];
    }

    double Options::GetDoubleVectorValue(std::string name, int nth, double def) {
        DoubleVector& vector = _doubleVectorMap[name];
        if (vector.size() > nth) {
            return vector[nth];
        }
        return def;
    }

    IntVector& Options::GetIntVector(std::string name) {
        return _intVectorMap[name];
    }

    int Options::GetIntVectorValue(std::string name, int nth, int def) {
        IntVector& vector = _intVectorMap[name];
        if (vector.size() > nth) {
            return vector[nth];
        }
        return def;
    }
    

    StringVector& Options::ParseOptions(int argc, char* argv[], CSimpleOpt::SOption* optionSpecs) {
        CSimpleOpt::SOption defaultSpecs[] = {
            SO_END_OF_OPTIONS
        };
        if (optionSpecs == NULL) {
            optionSpecs = defaultSpecs;
        }
        CSimpleOpt args(argc, argv, optionSpecs);
        while (args.Next()) {
            if (args.LastError() == SO_SUCCESS) {
                if (args.OptionArg()) {
                    this->SetString(args.OptionText(), string(args.OptionArg()));
                } else {
                    this->SetBool(args.OptionText(), true);
                }
            } else {
                // handle error (see the error codes - enum ESOError)
                std::cout << "ERROR:" << args.LastError() << std::endl;
            }
        }

        for (int i = 0; i < args.FileCount(); i++) {
            this->AppendString("args", args.File(i));
        }

        return this->GetStringVector("args");
    }

    std::ostream & operator<<(std::ostream &os, const Options& opt) {
        bool hasValue = false;
        for(Options::BoolMap::const_iterator iter = opt._boolMap.begin(); iter != opt._boolMap.end(); ++iter) {
            if (!hasValue) {
                os << "Options:";
                hasValue = true;
            }
            if (iter->first[0] != '_') {
                os << " " << (iter->second ? "+" : "-") << iter->first;
            }
        }
        if (hasValue) {
            os << endl;
            hasValue = false;
        }

        for(Options::IntMap::const_iterator iter = opt._intMap.begin(); iter != opt._intMap.end(); ++iter) {
            os << iter->first << " int " << iter->second << endl;
        }
        for(Options::DoubleMap::const_iterator iter = opt._doubleMap.begin(); iter != opt._doubleMap.end(); ++iter) {
            os << iter->first << " double " << iter->second << endl;
        }
        for(Options::StringMap::const_iterator iter = opt._stringMap.begin(); iter != opt._stringMap.end(); ++iter) {
            os << iter->first << " string " << iter->second << endl;
        }

        for(Options::IntVectorMap::const_iterator iter = opt._intVectorMap.begin(); iter != opt._intVectorMap.end(); ++iter) {
            IntVector vector = iter->second;
            os << iter->first << " ints";
            for (int i = 0; i < vector.size(); i++) {
                os << " " << vector[i];
            }
            os << endl;
        }

        for(Options::DoubleVectorMap::const_iterator iter = opt._doubleVectorMap.begin(); iter != opt._doubleVectorMap.end(); ++iter) {
            DoubleVector vector = iter->second;
            os << iter->first << " doubles";
            for (int i = 0; i < vector.size(); i++) {
                os << " " << vector[i];
            }
            os << endl;
        }

        for(Options::StringVectorMap::const_iterator iter = opt._stringVectorMap.begin(); iter != opt._stringVectorMap.end(); ++iter) {
            StringVector vector = iter->second;
            os << iter->first << " strings " << vector.size() << endl;
            for (int i = 0; i < vector.size(); i++) {
                os << vector[i] << endl;
            }
        }
        return os;
    }

    std::istream & operator>>(std::istream &is, Options& opt) {
        char *cbuf = new char[10240];
        while (is.good()) {
            is.getline(cbuf, 10240);
            if (is.good()) {
                stringstream ss(cbuf);
                string name;
                string type;
                ss >> name;
                // skip whitespace
                if (name[0] == '#' || name[0] == '%') {
                    continue;
                }
                if (name == "Options:") {
                    while (ss.good()) {
                        string option;
                        ss >> option;
                        bool value = (option[0] == '+');
                        string optionname = option.substr(1);
                        opt.SetBool(optionname, value);
                    }
                }
                ss >> type;
                if (name == OPTION_END) {
                    return is;
                }
                if (type == "int" ) {
                    int value;
                    ss >> value;
                    opt.SetInt(name, value);
                } else if (type == "double") {
                    double value;
                    ss >> value;
                    opt.SetDouble(name, value);
                } else if (type == "bool") {
                    string value;
                    ss >> value;
                    opt.SetBool(name, value == "true");
                } else if (type == "string") {
                    string value;
                    ss >> value;
                    opt.SetString(name, value);
                } else if (type == "ints") {
                    while (ss.good()) {
                        int val;
                        ss >> val;
                        opt.AppendInt(name, val);
                    }
                } else if (type == "doubles") {
                    while (ss.good()) {
                        double val;
                        ss >> val;
                        opt.AppendDouble(name, val);
                    }
                } else if (type == "strings") {
                    int len;
                    ss >> len;
                    for (int i = 0; i < len; i++) {
                        char sbuf[1024];
                        is.getline(sbuf, sizeof(sbuf));
                        if (is.good()) {
                            opt.AppendString(name, sbuf);
                        }
                    }
                }
            }
        }
        delete[] cbuf;
        return is;
    }
}

