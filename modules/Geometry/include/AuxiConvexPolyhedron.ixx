//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace merope {
template<unsigned short DIM>
inline vector<Segment<DIM>> edgesFromVertices(
    const vector<Point<DIM> >& renormalized_vertices,
    const vector<vector<long> >& face_indices) {
    std::set<array<long, 2>> indices_segment{};
    for (const auto& f_i : face_indices) {
        for (size_t j = 0; j + 1 < f_i.size(); j++) {
            indices_segment.insert({ min(f_i[j], f_i[j + 1]), max(f_i[j], f_i[j + 1]) });
        }
        indices_segment.insert({ min(f_i[0], f_i[f_i.size() - 1]), max(f_i[0], f_i[f_i.size() - 1]) });
    }
    //
    vector<Segment<DIM>> result{};
    for (const auto& i_seg : indices_segment) {
        result.emplace_back(Segment<DIM>{renormalized_vertices[i_seg[0]], renormalized_vertices[i_seg[1]]});
    }
    return result;
}


template<unsigned short DIM>
inline vector<HalfSpace<DIM> > facesFromVertices(
    const vector<Point<DIM>>& renormalized_vertices,
    const vector<vector<long>>& face_indices) {
    Point<DIM> insidePoint = average<DIM>(renormalized_vertices);
    vector<HalfSpace<DIM>> faces{};
    //
    for (const auto& face_index : face_indices) {
        Point<DIM> normal{};
        if constexpr (DIM == 3) {
            Point<DIM> vec1 = renormalized_vertices[face_index[1]] - renormalized_vertices[face_index[0]];
            Point<DIM> vec2 = renormalized_vertices[face_index[2]] - renormalized_vertices[face_index[1]];
            normal = geomTools::prodVec<DIM>(vec1, vec2);
        } else {
            Point<DIM> vec1 = renormalized_vertices[face_index[1]] - renormalized_vertices[face_index[0]];
            normal = { vec1[1], -vec1[0] };
        }
        geomTools::renormalize<DIM>(normal);
        HalfSpace<DIM> hspace(normal, renormalized_vertices[face_index[0]]);
        if (not hspace.isInside(insidePoint)) {
            hspace = HalfSpace<DIM>(-1 * normal, renormalized_vertices[face_index[0]]);
        }
        faces.push_back(hspace);
    }
    return faces;
}
}  // namespace merope




