//
//  piOptionParser.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/24/13.
//
//

#include <fstream>
#include <sstream>
#include "piOptionParser.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filestream.h"

namespace pi {
#pragma mark Private Members

    std::string OptionParser::stringValue(std::string name) {
        if (_jsonDoc.IsArray()) {
            for (int i = 0; i < _jsonDoc.Size(); i++) {
                if (_jsonDoc[i].IsObject()) {
                    if (_jsonDoc[i].HasMember(name.c_str())) {
                        return stringValue(_jsonDoc[i][name.c_str()]);
                    }
                }
            }
        } else {
            return stringValue(_jsonDoc[name.c_str()]);
        }
        return "";
    }

#pragma mark Public Members
    std::string OptionParser::json(std::vector<int>& data) {
        std::stringstream out;
        out << "[";
        for (int i = 0; i < data.size(); i++) {
            if (i > 0) {
                out << ", ";
            }
            out << data[i];
        }
        out << "]";
        return out.str();
    }

    std::string OptionParser::json(std::vector<double>& data) {
        std::stringstream out;
        out << "[";
        for (int i = 0; i < data.size(); i++) {
            if (i > 0) {
                out << ", ";
            }
            out << data[i];
        }
        out << "]";
        return out.str();
    }

    std::string OptionParser::json(std::vector<std::string>& data) {
        std::stringstream out;
        out << "[";
        for (int i = 0; i < data.size(); i++) {
            if (i > 0) {
                out << ", ";
            }
            out << data[i];
        }
        out << "]";
        return out.str();
    }


    void OptionParser::value(std::string& name, int& var) {
        var = _jsonDoc[name.c_str()].GetInt();
    }

    void OptionParser::value(std::string& name, double& var) {
        var = _jsonDoc[name.c_str()].GetDouble();
    }

    void OptionParser::value(std::string& name, std::string& var) {
        var = _jsonDoc[name.c_str()].GetString();
    }


    void OptionParser::value(int idx, unsigned long& var) {
        var = _jsonDoc[idx].GetInt();
    }

    void OptionParser::value(int idx, unsigned int& var) {
        var = _jsonDoc[idx].GetInt();
    }

    void OptionParser::value(int idx, signed long& var) {
        var = _jsonDoc[idx].GetInt();
    }

    void OptionParser::value(int idx, int& var) {
        var = _jsonDoc[idx].GetInt();
    }

    void OptionParser::value(int idx, double& var) {
        var = _jsonDoc[idx].GetDouble();
    }

    void OptionParser::value(int idx, std::string& var) {
        var = _jsonDoc[idx].GetString();
    }


    void OptionParser::values(std::string& name, std::vector<int>& var) {
        values(_jsonDoc[name.c_str()], var);
    }

    void OptionParser::values(std::string& name, std::vector<double>& var) {
        values(_jsonDoc[name.c_str()], var);
    }

    void OptionParser::values(std::string& name, std::vector<std::string>& var) {
        values(_jsonDoc[name.c_str()], var);
    }

    void OptionParser::values(int idx, std::vector<int>& var) {
        values(_jsonDoc[idx], var);
    }

    void OptionParser::values(int idx, std::vector<double>& var) {
        values(_jsonDoc[idx], var);
    }

    void OptionParser::values(int idx, std::vector<std::string>& var) {
        values(_jsonDoc[idx], var);
    }

    void OptionParser::values(rapidjson::Value& arrayValue, std::vector<int>& var) {
        if (!arrayValue.IsArray()) {
            return;
        }
        for (int i = 0; i < arrayValue.Size(); i++) {
            var.push_back(arrayValue[i].GetInt());
        }

    }

    void OptionParser::values(rapidjson::Value& arrayValue, std::vector<double>& var) {
        if (!arrayValue.IsArray()) {
            return;
        }
        for (int i = 0; i < arrayValue.Size(); i++) {
            var.push_back(arrayValue[i].GetDouble());
        }
    }

    void OptionParser::values(rapidjson::Value& arrayValue, std::vector<std::string>& var) {
        if (!arrayValue.IsArray()) {
            return;
        }
        for (int i = 0; i < arrayValue.Size(); i++) {
            var.push_back(arrayValue[i].GetString());
        }
    }

    std::string OptionParser::stringValue() {
        return stringValue(_jsonDoc);
    }

    std::string OptionParser::stringValue(const rapidjson::Document::ValueType& value) {
        std::stringstream out;
        if (value.IsNumber()) {
            out << value.GetDouble();
        } else if (value.IsString()) {
            out << value.GetString();
        } else if (value.IsArray()) {
            out << stringValues(value);
        } else if (value.IsObject()) {
            rapidjson::Document::ConstMemberIterator iter = value.MemberBegin();
            while (iter != value.MemberEnd()) {
                out << iter->name.GetString() << ":" << stringValue(iter->value);
                iter++;
            }
        }
        return out.str();
    }


    std::string OptionParser::stringValues(const rapidjson::Document::ValueType& value) {
        std::stringstream out;
        if (value.IsObject()) {
        } else if (value.IsArray()) {
            out << "[";
            for (int i = 0; i < value.Size(); i++) {
                if (i > 0) {
                    out << ", ";
                }
                out << stringValue(value[i]);
            }
            out << "]";
        } else {
            out << stringValue(value);
        }
        return out.str();
    }

    void OptionParser::read(std::string filename) {
        FILE* f = fopen(filename.c_str(), "r");
        rapidjson::FileStream fs(f);
        _jsonDoc.ParseStream<0>(fs);
        if (_jsonDoc.IsArray()) {
            _size = _jsonDoc.Size();
        } else {
            _size = 1;
        }
        fclose(f);
    }

    void OptionParser::write(std::string filename) {
        using namespace rapidjson;
        FileStream f(stdout);
        PrettyWriter<FileStream> writer(f);
        _jsonDoc.Accept(writer);
    }
}