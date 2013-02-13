//
//  myOptions.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "piOptions.h"

namespace pi {
    void Options::Set(std::string name, bool value) {
        std::pair<BoolMap::iterator, bool> result = _boolMap.insert(BoolPair(name, value));
        if (!result.second) {
            result.first->second = value;
        }
    }
    
    void Options::Set(std::string name, int value) {
        std::pair<IntMap::iterator, bool> result = _intMap.insert(IntPair(name, value));
        if (!result.second) {
            result.first->second = value;
        }
    }

    void Options::Set(std::string name, double value) {
        std::pair<DoubleMap::iterator, bool> result = _doubleMap.insert(DoublePair(name, value));
        if (!result.second) {
            result.first->second = value;
        }    }

    void Options::Set(std::string name, std::string value) {
        std::pair<StringMap::iterator, bool> result = _stringMap.insert(StringPair(name, value));
        if (!result.second) {
            result.first->second = value;
        }
    }

    void Options::Add(std::string name, std::string value) {
        _stringVectorMap[name].push_back(value);
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
                    this->Set(args.OptionText(), args.OptionArg());
                } else {
                    this->Set(args.OptionText(), true);
                }
            } else {
                // handle error (see the error codes - enum ESOError)
                std::cout << "ERROR:" << args.LastError() << std::endl;
            }
        }

        for (int i = 0; i < args.FileCount(); i++) {
            this->Add("args", args.File(i));
        }
        return this->GetStringVector("args");
    }

}