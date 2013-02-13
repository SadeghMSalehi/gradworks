//
//  myOptions.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#ifndef __ParticlesGUI__myOptions__
#define __ParticlesGUI__myOptions__

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "SimpleOpt.h"

namespace pi {
    typedef std::vector<std::string> StringVector;
    
    class Options {
    private:
        typedef std::map<std::string, std::string> StringMap;
        typedef std::pair<std::string, std::string> StringPair;
        typedef std::map<std::string, int> IntMap;
        typedef std::pair<std::string, int> IntPair;
        typedef std::map<std::string, double> DoubleMap;
        typedef std::pair<std::string, double> DoublePair;
        typedef std::map<std::string, bool> BoolMap;
        typedef std::pair<std::string, bool> BoolPair;
        typedef std::map<std::string, StringVector> StringVectorMap;
        typedef std::pair<std::string, StringVector> StringVectorPair;
        
    public:
        void Set(std::string name, bool value);
        void Set(std::string name, int value);
        void Set(std::string name, double value);
        void Set(std::string name, std::string value);
        void Add(std::string name, std::string value);

        bool GetBool(std::string name, bool def = false);
        int GetInt(std::string name, int def);
        int GetStringAsInt(std::string name, int def);
        double GetDouble(std::string name, double def);
        double GetStringAsDouble(std::string name, double def);
        std::string GetString(std::string name, std::string def);
        StringVector& GetStringVector(std::string name);
        StringVector& ParseOptions(int argc, char* argv[], CSimpleOpt::SOption*);

    private:
        BoolMap _boolMap;
        StringMap _stringMap;
        IntMap _intMap;
        DoubleMap _doubleMap;
        StringVectorMap _stringVectorMap;
    };
}
#endif /* defined(__ParticlesGUI__myOptions__) */
