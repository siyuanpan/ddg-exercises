// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:

    for (size_t i = 0; i < mesh->nVertices(); i++) {
        geometry->vertexIndices[i] = i;
    }

    // This is also a valid way to access indices (directly by mesh element).
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        geometry->vertexIndices[v] = idx;
        idx++;
    }

    for (size_t i = 0; i < mesh->nEdges(); i++) {
        geometry->edgeIndices[i] = i;
    }

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        geometry->faceIndices[i] = i;
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)

    // You can more easily get the indices of mesh elements using the function getIndex(), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh->nEdges()*2);
    for (Edge e : mesh->edges()) {
        tripletList.push_back(T(geometry->edgeIndices[e], geometry->vertexIndices[e.firstVertex()], 1));
        tripletList.push_back(T(geometry->edgeIndices[e], geometry->vertexIndices[e.secondVertex()], 1));
    }
    SparseMatrix<size_t> mat(mesh->nEdges(), mesh->nVertices());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat; // placeholder
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh->nFaces()*3);
    for (Face f : mesh->faces()) {
        for (Edge e : f.adjacentEdges()) {
            tripletList.push_back(T(geometry->faceIndices[f], geometry->edgeIndices[e], 1));
        }
    }
    SparseMatrix<size_t> mat(mesh->nFaces(), mesh->nEdges());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat; // placeholder
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nVertices());
    for (auto i : subset.vertices) {
        vec[i] = 1;
    }
    return vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nEdges());
    for (auto i : subset.edges) {
        vec[i] = 1;
    }
    return vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nFaces());
    for (auto i : subset.faces) {
        vec[i] = 1;
    }
    return vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    vertices.insert(subset.vertices.begin(), subset.vertices.end());
    edges.insert(subset.edges.begin(), subset.edges.end());
    faces.insert(subset.faces.begin(), subset.faces.end());

    auto A0 = buildVertexEdgeAdjacencyMatrix();
    auto A1 = buildFaceEdgeAdjacencyMatrix();
    auto vertexVec = buildVertexVector(subset);
    auto edgeVec = buildEdgeVector(subset);

    auto starVertexEdge = (A0 * vertexVec).eval();
    for (int i = 0; i < starVertexEdge.size(); ++i) {
        if (starVertexEdge(i)) {
            edges.insert(i);
        }
    }

    auto starVertexFace = (A1 * A0 * vertexVec).eval();
    for (int i = 0; i < starVertexFace.size(); ++i) {
        if (starVertexFace(i)) {
            faces.insert(i);
        }
    }

    auto startEdgeFace = (A1 * edgeVec).eval();
    for (int i = 0; i < startEdgeFace.size(); ++i) {
        if (startEdgeFace(i)) {
            faces.insert(i);
        }
    }
  
    MeshSubset S = MeshSubset(vertices, edges, faces);
    return S; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    vertices.insert(subset.vertices.begin(), subset.vertices.end());
    edges.insert(subset.edges.begin(), subset.edges.end());
    faces.insert(subset.faces.begin(), subset.faces.end());

    auto A0 = buildVertexEdgeAdjacencyMatrix();
    auto A1 = buildFaceEdgeAdjacencyMatrix();
    auto edgeVec = buildEdgeVector(subset);
    auto faceVec = buildFaceVector(subset);

    auto starFaceEdge = (A1.transpose() * faceVec).eval();
    for (int i = 0; i < starFaceEdge.size(); ++i) {
        if (starFaceEdge(i)) {
            edges.insert(i);
        }
    }

    auto starFaceVertex = ((A1*A0).transpose() * faceVec).eval();
    for (int i = 0; i < starFaceVertex.size(); ++i) {
        if (starFaceVertex(i)) {
            vertices.insert(i);
        }
    }

    auto startEdgeVertex = (A0.transpose() * edgeVec).eval();
    for (int i = 0; i < startEdgeVertex.size(); ++i) {
        if (startEdgeVertex(i)) {
            vertices.insert(i);
        }
    }
  
    MeshSubset S = MeshSubset(vertices, edges, faces);
    return S; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    MeshSubset S = closure(star(subset));
    MeshSubset C = star(closure(subset));

    S.deleteSubset(C);
    
    return S; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    for (auto i : subset.faces) {
        Face f = mesh->face(i);
        for (Edge e : f.adjacentEdges()) {
            if (subset.edges.find(geometry->edgeIndices[e]) == subset.edges.end()) return false;
        }
        for (Vertex v : f.adjacentVertices()) {
            if (subset.vertices.find(geometry->vertexIndices[v]) == subset.vertices.end()) return false;
        }
    }
    for (auto i : subset.edges) {
        Edge e = mesh->edge(i);
        for (Vertex v : e.adjacentVertices()) {
            if (subset.vertices.find(geometry->vertexIndices[v]) == subset.vertices.end()) return false;
        }
    }
    return true; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    MeshSubset S = subset.deepCopy();
    std::set<size_t> delEdges;
    std::set<size_t> delVertices;
    for (auto i : S.faces) {
        Face f = mesh->face(i);
        for (Edge e : f.adjacentEdges()) {
            if (S.edges.find(geometry->edgeIndices[e]) == S.edges.end()) return -1;
            delEdges.insert(geometry->edgeIndices[e]);
        }
        for (Vertex v : f.adjacentVertices()) {
            if (S.vertices.find(geometry->vertexIndices[v]) == S.vertices.end()) return -1;
            delVertices.insert(geometry->vertexIndices[v]);
        }
    }
    S.deleteEdges(delEdges);
    for (auto i : S.edges) {
        Edge e = mesh->edge(i);
        for (Vertex v : e.adjacentVertices()) {
            if (S.vertices.find(geometry->vertexIndices[v]) == S.vertices.end()) return -1;
            delVertices.insert(geometry->vertexIndices[v]);
        }
    }
    S.deleteVertices(delVertices);
    if (!S.faces.empty() && S.edges.empty() && S.vertices.empty()) return 2;
    else if (S.faces.empty() && !S.edges.empty() && S.vertices.empty()) return 1;
    else if (S.faces.empty() && S.edges.empty() && !S.vertices.empty()) return 0;
    else return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    std::set<size_t> vertices;
    std::set<size_t> edges;
    int degree = isPureComplex(subset);
    assert(degree != -1);
    if (degree == 2) {
        for (auto i : subset.edges) {
            int cnt = 0;
            for (Face f : mesh->edge(i).adjacentFaces()) {
                if (subset.faces.find(geometry->faceIndices[f]) != subset.faces.end()) cnt++;
            }
            if (cnt != 2) {
                edges.insert(i);
                vertices.insert(geometry->vertexIndices[mesh->edge(i).firstVertex()]);
                vertices.insert(geometry->vertexIndices[mesh->edge(i).secondVertex()]);
            }
        }
    } else if (degree == 1) {
        for (auto i : subset.vertices) {
            int cnt = 0;
            for (Edge e : mesh->vertex(i).adjacentEdges()) {
                if (subset.edges.find(geometry->edgeIndices[e]) != subset.edges.end()) cnt++;
            }
            if (cnt != 2) vertices.insert(i);
        }
    }

    MeshSubset S = MeshSubset(vertices, edges, {});
    return S; // placeholder
}