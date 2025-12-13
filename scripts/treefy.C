// see.C
#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <TString.h>
#include <TRegexp.h>
#include <iostream>
#include <fstream>
#include <set>
#include <regex>

void treefy(const char* infile = "test.root",
            std::string outfile = "",
            int nevt = 0) {
    std::string infile_str(infile);
    if (outfile == "") {
      outfile = infile_str.substr(0, infile_str.find_last_of(".")) + "_ttree.root";
    }

    ROOT::RDataFrame df("Events", infile);

    // Collect column names and types
    auto colNames = df.GetColumnNames();
    std::set<std::string> classes;
    std::regex class_pattern(R"(ROOT::[^\n]+)");

    for (auto&& col : colNames) {
        std::string colType = df.GetColumnType(col);
        std::smatch match;
        if (std::regex_search(colType, match, class_pattern)) {
            std::string cls = match.str(0);
            if (std::count(cls.begin(), cls.end(), '<') != std::count(cls.begin(), cls.end(), '>')) {
                while (!cls.empty() && cls.back() == '>') cls.pop_back();
            }
            classes.insert(cls);
        }
    }

    if (gSystem->AccessPathName("dict.so")) {
      std::cout << "Building dict" << std::endl;
      // Write dict.h
      std::ofstream fout("dict.h");
      fout << "#include \"ROOT/RVec.hxx\"\n";
      fout << "#include \"Math/Vector4D.h\"\n";
      fout << "#include \"Math/GenVector/DisplacementVector3D.h\"\n";
      fout << "#include \"Math/GenVector/Cartesian3D.h\"\n";
      fout << "#include \"Math/GenVector/PositionVector3D.h\"\n";
      fout << "#include \"Math/GenVector/VectorUtil.h\"\n";
      fout << "#include \"Math/SVector.h\"\n";
      fout << "#include \"Math/SMatrix.h\"\n";
      fout << "#ifdef __CLING__\n";
      for (const auto& cls : classes) {
        fout << "#pragma link C++ class " << cls << "+;\n";
      }
      fout << "#endif\n";
      fout.close();

      // Build dictionary
      gSystem->Exec("rootcling -f dict.cxx -c dict.h");
      gSystem->Exec("g++ -shared -fPIC $(root-config --cflags --libs) dict.cxx -o dict.so");
    }
    gSystem->Load("dict.so");

    // Select columns
    std::vector<std::string> filtered;
    // std::regex part_pattern("^Part_fourMomentum");
    std::vector<std::string> bad_types = {"Event_shortDstVersion"};

    for (auto&& col : colNames) {
        std::string colstr = std::string(col);
        bool bad = false;
        for (auto&& badpat : bad_types) {
          if (colstr.find(badpat) != std::string::npos) {
            bad = true;
            break;
          }
        }
        if (!bad) filtered.push_back(colstr);
    }


    // Snapshot
    if (nevt > 0) {
      df.Range(nevt).Snapshot("t", outfile, filtered);
    } else {
      df.Snapshot("t", outfile, filtered);
    }
}
