// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplet;
    triplet.reserve(mesh.nVertices());
    for (Vertex v : mesh.vertices()) {
        triplet.push_back(T(vertexIndices[v], vertexIndices[v], barycentricDualArea(v)));
    }
    SparseMatrix<double> mat(mesh.nVertices(), mesh.nVertices());
    mat.setFromTriplets(triplet.begin(), triplet.end());
    return mat; // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplet;
    triplet.reserve(mesh.nEdges());
    for (Edge e : mesh.edges()) {
        Halfedge h1 = e.halfedge();
        Halfedge h2 = h1.twin();
        triplet.push_back(T(edgeIndices[e], edgeIndices[e], (cotan(h1)+cotan(h2))*0.5));
    }
    SparseMatrix<double> mat(mesh.nEdges(), mesh.nEdges());
    mat.setFromTriplets(triplet.begin(), triplet.end());
    return mat; // placeholder
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplet;
    triplet.reserve(mesh.nFaces());
    for (Face f : mesh.faces()) {
        triplet.push_back(T(faceIndices[f], faceIndices[f], 1./ faceArea(f)));
    }
    SparseMatrix<double> mat(mesh.nFaces(), mesh.nFaces());
    mat.setFromTriplets(triplet.begin(), triplet.end());
    return mat; // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplet;
    triplet.reserve(mesh.nEdges()*2);
    for (Edge e : mesh.edges()) {
        Halfedge hf = e.halfedge();
        triplet.push_back(T(edgeIndices[e], vertexIndices[hf.tipVertex()], 1.));
        triplet.push_back(T(edgeIndices[e], vertexIndices[hf.tailVertex()], -1.));
    }
    SparseMatrix<double> mat(mesh.nEdges(), mesh.nVertices());
    mat.setFromTriplets(triplet.begin(), triplet.end());
    return mat; // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplet;
    triplet.reserve(mesh.nFaces()*3);
    for (Face f : mesh.faces()) {
        for (Halfedge hf : f.adjacentHalfedges()) {
            Edge e = hf.edge();
            Halfedge h = e.halfedge();
            triplet.push_back(T(faceIndices[f], edgeIndices[hf.edge()], (hf == h)?1.:-1.));
        }
    }
    SparseMatrix<double> mat(mesh.nFaces(), mesh.nEdges());
    mat.setFromTriplets(triplet.begin(), triplet.end());
    return mat; // placeholder
}

} // namespace surface
} // namespace geometrycentral