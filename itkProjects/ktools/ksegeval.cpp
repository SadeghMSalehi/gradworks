//
//  ksegeval.cpp
//  ktools
//
//  Created by Joohwi Lee on 3/7/13.
//
//

#include "ksegeval.h"

#include "piImageDef.h"
#include "piImageIO.h"
#include "piOptions.h"
#include "piImageProcessing.h"


using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    __noverbose = 1;
    CSimpleOpt::SOption specs[] = {
        { 0, "-o", SO_REQ_SEP },
        { 1, "--template", SO_REQ_SEP },
        { 2, "--columns", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };

    Options argParser;
    StringVector args = argParser.ParseOptions(argc, argv, specs);    string eq = argParser.GetString("-e");
    string outputFilename = argParser.GetString("-o");
    string templ = argParser.GetString("--template");

    if (args.size() == 0) {
        cout << argv[0] << " -o output-csv [--template template segmentation] [ --columns A,B,AB,dice,total] input1 input2 ... (standard1 standard2 ...)" << endl;
        cout << "          if --template is not given, the arguments must be in a pair (ex. input1 input2 standard1 standard2) " << endl;
        cout << "          if --template is given, every argument is compared to the template." << endl;
        cout << "          --columns will accept the order of the columns (A,B,A and B, total, dice)" << endl;
        return 0;
    }

    ImageIO<LabelImage> io;
    std::vector<AtlasSimilarityScore> scores;
    if (templ == "") {
        int nsubj = args.size() / 2;
        for (int i = 0; i < nsubj; i++) {
            AtlasSimilarityScore score;
            score.Compute(io.ReadCastedImage(args[i]), io.ReadCastedImage(args[i+nsubj]));
            scores.push_back(score);
        }
    } else {
        LabelImage::Pointer tmplImg = io.ReadCastedImage(templ);
        for (int i = 0; i < args.size(); i++) {
            AtlasSimilarityScore score;
            score.Compute(io.ReadCastedImage(args[i]), tmplImg);
            scores.push_back(score);
        }
    }

    ostream* os = &cout;
    ofstream fout;
    if (outputFilename != "") {
        fout.open(outputFilename.c_str());
        os = &fout;
    }
    StringVector cols = argParser.GetSplitString("--columns", ",", "A,AB,B,dice");
    for (int j = 0; j < scores[0].labelMap.size(); j++) {
        for (int i = 0; i < scores.size(); i++) {
            for (int k = 0; k < cols.size(); k++) {
                if (cols[k] == "A") (*os) << scores[i].labelMap[j].A;
                if (cols[k] == "B") (*os) << scores[i].labelMap[j].B;
                if (cols[k] == "AB") (*os) << scores[i].labelMap[j].AB;
                if (cols[k] == "dice") (*os) << scores[i].labelMap[j].dice();
                if (cols[k] == "total") (*os) << scores[i].labelMap[j].total();
                (*os) << " ";
            }
        }
        (*os) << endl;
    }
    if (fout.is_open()) {
        fout.close();
    }
    os = NULL;
}