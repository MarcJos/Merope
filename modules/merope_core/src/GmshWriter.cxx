//! Copyright : see license.txt
//!
//! \brief
//


#include "Mesh/GmshWriter.hxx"


namespace merope {

void mesh::gmsh_writer::auxi::writeEnd(std::ostream& f, vector<string> name_outputs, double meshSize, size_t meshOrder, bool binaryOutput) {
    f << endl;
    f << "MeshSize {:} =" << meshSize << ";\n";
    f << endl;
    f << "Mesh.ElementOrder =" << meshOrder << ";\n";
    f << "Mesh 2;\n";
    f << "Mesh 3;\n";
    if (binaryOutput) {
        f << "Mesh.Binary=1;\n";
    }
    // save
    f << "Mesh.MshFileVersion = 2.2;\n";
    for (auto name : name_outputs) {
        f << R"(Save ")" << name << R"(";)" << "\n";
    }
}

}  // namespace merope

