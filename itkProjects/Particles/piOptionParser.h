//
//  piOptionParser.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/24/13.
//
//

#ifndef __ParticleGuidedRegistration__piOptionParser__
#define __ParticleGuidedRegistration__piOptionParser__

#include <iostream>
#include <string>
#include <vector>

#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"


namespace pi {
    class OptionParser {

    public:
        OptionParser() {}
        OptionParser(int argc, const char* argv[]) {}
        OptionParser(std::string json) {
            _jsonDoc.Parse<0>(json.c_str());
            if (_jsonDoc.IsArray()) {
                _size = _jsonDoc.Size();
            } else {
                _size = 1;
            }
        }

        inline const int size() { return _size; }
        inline const int size(std::string name) {
            return _jsonDoc[name.c_str()].Size();
        }
        inline const int size(int idx) {
            return _jsonDoc[idx].Size();
        }

        void read(std::string filename);
        void write(std::string filename);

        std::string stringValue();
        std::string stringValue(std::string name);

        std::string json(std::vector<int>& data);
        std::string json(std::vector<double>& data);
        std::string json(std::vector<std::string>& data);

        void value(std::string& name, int& var);
        void value(std::string& name, double& var);
        void value(std::string& name, std::string& var);

        void value(int idx, unsigned int& var);
        void value(int idx, unsigned long& var);
        void value(int idx, signed long& var);

        void value(int idx, int& var);
        void value(int idx, double& var);
        void value(int idx, std::string& var);

        void values(std::string& name, std::vector<int>& var);
        void values(std::string& name, std::vector<double>& var);
        void values(std::string& name, std::vector<std::string>& var);

        void values(int idx, std::vector<int>& var);
        void values(int idx, std::vector<double>& var);
        void values(int idx, std::vector<std::string>& var);

        template<typename S>
        void values(S* vars, int nVars, int nElems) {
            int k = 0;
            for (int i = 0; i < nVars; i++) {
                for (int j = 0; j < nElems; j++) {
                    double d = 0;
                    value(k, d);
                    vars[i][j] = d;
                    k++;
                }
            }
        }

    private:
        void values(rapidjson::Value& value, std::vector<int>& var);
        void values(rapidjson::Value& value, std::vector<double>& var);
        void values(rapidjson::Value& value, std::vector<std::string>& var);


        std::string stringValue(const rapidjson::Document::ValueType& value);
        std::string stringValues(const rapidjson::Document::ValueType& value);

        int _size;
        rapidjson::Document _jsonDoc;
    };


}
#endif /* defined(__ParticleGuidedRegistration__piOptionParser__) */