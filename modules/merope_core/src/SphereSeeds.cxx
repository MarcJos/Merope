//! Copyright : see license.txt
//!
//! \brief

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "MultiInclusions/SphereSeeds.hxx"


#include "MeropeNamespace.hxx"


namespace merope {


void writeLongLine(const vector<string>& l, ostream& fout) {
    string lp = "", ls = "";
    for (auto& i : l) {
        if (lp.size() == 0 and ls.size() > 0) lp = ls;
        if ((lp + i).size() > 72) {
            fout << lp << endl;
            lp = "";
            ls = i;
        } else {
            lp += i;
            ls = "";
        }
    }
    if (lp.size() > 0) fout << lp << endl;
    if (ls.size() > 0) fout << ls << endl;
}

void readKeyword(const std::string keyword, ifstream& f) {
    string str;
    f >> str;
    if (keyword != str) {
        throw logic_error(
            "Read [" + str + "] while waiting for keyword [" + keyword
            + "]");
    }
}

} // namespace merope

