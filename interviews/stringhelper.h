//
//  stringhelper.h
//  interviews
//
//  Created by Joowhi Lee on 5/24/15.
//
//

#ifndef interviews_stringhelper_h
#define interviews_stringhelper_h

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>

struct StringHelper {
    
    /*
     This function splits a string into substrings with respect to the given token. The token could be multiple letters. If noempty is true, any empty string produced by the split is ignored in the output.
     */
    std::vector<std::string> split(std::string txt, std::string tok, bool noempty = false) {
        std::vector<std::string> list;
        for (int j = 0; j < txt.size(); ) {
            size_t n = txt.find(tok, j);
            if (n == std::string::npos) {
                if (!(noempty && n-j == 0)) {
                    list.push_back(txt.substr(j, n-j));
                }
                break;
            } else {
                if (!(noempty && n-j == 0)) {
                    list.push_back(txt.substr(j, n-j));
                }
                j = n + tok.size();
            }
        }
        return list;
    }
	
	std::vector<int> split2int(std::string txt, std::string tok, bool noempty = false) {
		return convert2int(split(txt, tok, noempty));
	}
	
    
    std::vector<std::vector<std::string> > lstlst(std::string s, std::string tok1, std::string tok2) {
		std::vector<std::vector<std::string> > itemlist;
		std::vector<std::string> seglist = split(s, tok1, true);
		for (int j = 0; j < seglist.size(); j++) {
			itemlist.push_back(split(seglist[j], tok2, true));
		}
		return itemlist;
    }
	
	std::vector<std::vector<int> > lstlst2int(std::string s, std::string tok1, std::string tok2) {
		std::vector<std::vector<int> > itemlist;
		std::vector<std::string> seglist = split(s, tok1, true);
		for (int j = 0; j < seglist.size(); j++) {
			itemlist.push_back(convert2int(split(seglist[j], tok2, true)));
		}
		return itemlist;
	}
	
	std::vector<std::vector<double> > lstlst2dbl(std::string s, std::string tok1, std::string tok2) {
		std::vector<std::vector<double> > itemlist;
		std::vector<std::string> seglist = split(s, tok1, true);
		for (int j = 0; j < seglist.size(); j++) {
			itemlist.push_back(convert2dbl(split(seglist[j], tok2, true)));
		}
		return itemlist;
	}
	
    std::vector<double> convert2dbl(std::vector<std::string> lst) {
        std::vector<double> e;
        for (int j = 0; j < lst.size(); j++) {
            e.push_back(str2dbl(lst[j]));
        }
        return e;
    }

    std::vector<int> convert2int(std::vector<std::string> lst) {
        std::vector<int> e;
        for (int j = 0; j < lst.size(); j++) {
            e.push_back(str2dbl(lst[j]));
        }
        return e;
    }

    
    int str2int(std::string s) {
        return atoi(s.c_str());
    }
    
    double str2dbl(std::string s) {
        return atof(s.c_str());
    }
	
	template<typename T>
    void print(std::vector<T>& lst) {
        for (int j = 0; j < lst.size(); j++) {
            std::cout << lst[j] << ",";
        }
        std::cout << std::endl;
    }

	template<typename T>
	void print(std::vector<std::vector<T> >& lst) {
		for (int j = 0; j < lst.size(); j++) {
			print(lst[j]);
		}
		std::cout << std::endl;
	}
};

#endif
