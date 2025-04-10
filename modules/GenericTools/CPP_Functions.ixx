//! Copyright : see license.txt
//!
//! \briefPurely informatic function for doing some C++ manipulations
//

#pragma once


namespace merope {
namespace cppFunctions {

template<typename T>
class UnionFind {
public:
    unordered_map<T, T> parent;

    void add(const T& x) {
        if (parent.find(x) == parent.end()) {
            parent[x] = x;
        }
    }

    T find(const T& x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    void unite(const T& x, const T& y) {
        T rootX = find(x);
        T rootY = find(y);
        if (rootX != rootY) {
            parent[rootX] = rootY;
        }
    }
};

template<typename T>
vector<unordered_set<T>> computeEquivalenceClass(const vector<pair<T, T>>& pairs) {
    UnionFind<T> uf;
    for (const auto& p : pairs) {
        uf.add(p.first);
        uf.add(p.second);
        uf.unite(p.first, p.second);
    }

    unordered_map<T, int> rootToClass;
    vector<unordered_set<T>> classes;
    for (const auto& p : pairs) {
        T root = uf.find(p.first);
        if (rootToClass.find(root) == rootToClass.end()) {
            unordered_set<T> newSet;
            classes.push_back(newSet);
            rootToClass[root] = classes.size() - 1;
        }
        classes[rootToClass[root]].insert(p.first);
        classes[rootToClass[root]].insert(p.second);
    }

    return classes;
}


}  // namespace cppFunctions
}  // namespace merope