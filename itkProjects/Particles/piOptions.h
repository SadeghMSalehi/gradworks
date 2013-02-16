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

/**
 * Todo
 * 
 *   value interpolated string 
 *      example)
 *         NumberOfParticles: int 100
 *         Name: interpolated_string file_#NumberOfParticles#.txt
 *         GetString("Name", "") => file_100.txt
 */
#define OPTION_END "==option_end=="
namespace pi {
    typedef std::vector<std::string> StringVector;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<int> IntVector;
    
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
        typedef std::map<std::string, DoubleVector> DoubleVectorMap;
        typedef std::pair<std::string, DoubleVector> DoubleVectorPair;
        typedef std::map<std::string, IntVector> IntVectorMap;
        typedef std::pair<std::string, IntVector> IntVectorPair;
    public:
        void Set(std::string name, bool value);
        void Set(std::string name, int value);
        void Set(std::string name, double value);
        void Set(std::string name, std::string value);
        void AppendString(std::string name, std::string value);
        void AppendDouble(std::string name, double value);
        void AppendInt(std::string name, int value);

        bool GetBool(std::string name, bool def = false);
        int GetInt(std::string name, int def);
        int GetStringAsInt(std::string name, int def);
        double GetDouble(std::string name, double def);
        double GetStringAsDouble(std::string name, double def);
        std::string GetString(std::string name, std::string def);

        StringVector& GetStringVector(std::string name);
        std::string GetStringVectorValue(std::string name, int i, std::string def = "");
        
        DoubleVector& GetDoubleVector(std::string name);
        double GetDoubleVectorValue(std::string name, int nth, double def = 0);

        IntVector& GetIntVector(std::string name);
        int GetIntVectorValue(std::string name, int nth, int def = 0);

        StringVector& ParseOptions(int argc, char* argv[], CSimpleOpt::SOption*);

    private:
        BoolMap _boolMap;
        StringMap _stringMap;
        IntMap _intMap;
        DoubleMap _doubleMap;
        StringVectorMap _stringVectorMap;
        DoubleVectorMap _doubleVectorMap;
        IntVectorMap _intVectorMap;

        friend std::ostream & operator<<(std::ostream &os, const Options& opt);
        friend std::istream & operator>>(std::istream &is, Options& opt);
    };
}
#endif /* defined(__ParticlesGUI__myOptions__) */
