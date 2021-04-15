// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge.
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

    // TODO
    if (he.isInterior()) {
        Halfedge h1 = he.next().next();
        Halfedge h2 = he.next().twin();
        Vector3 v1 = halfedgeVector(h1);
        Vector3 v2 = halfedgeVector(h2);
        return dot(v1, v2)/norm(cross(v1, v2));
    } else {
        return 0; // placeholder
    }
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    // TODO
    double inv3 = 1./3.;
    double dualArea = 0;
    if (v.isBoundary()) {
        for (Halfedge he : v.outgoingHalfedges()) {
            if (he.isInterior()) {
                dualArea += faceArea(he.face())*inv3;
            }
        }
    } else {
        for (Face f : v.adjacentFaces()) {
            dualArea += faceArea(f)*inv3;
        }
    }
    return dualArea; // placeholder
}

/*
 * Computes the angle (in radians) at a given corner.
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    // TODO
    Halfedge he = c.halfedge();
    Vector3 v1 = halfedgeVector(he);
    he = he.next().next().twin();
    Vector3 v2 = halfedgeVector(he);
    double q = dot(unit(v1), unit(v2));
    q = clamp(q, -1., 1.);
    double angle = std::acos(q);
    return angle; // placeholder
}

/*
 * Computes the signed angle (in radians) between two adjacent faces.
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    // TODO
    Vector3 v = halfedgeVector(he);
    Vector3 n1 = faceNormal(he.face());
    Vector3 n2 = faceNormal(he.twin().face());
    double theta = std::atan2(dot(unit(v), cross(n1, n2)), dot(n1, n2));
    return theta; // placeholder
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    // TODO
    Vector3 normal{0, 0, 0};
    for (Face f : v.adjacentFaces()) {
        Vector3 n = faceNormal(f);
        normal += n;
    }
    return unit(normal); // placeholder
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    // TODO
    Vector3 normal{0, 0, 0};
    for (Corner c : v.adjacentCorners()) {
        Vector3 n = faceNormal(c.face());
        double phi = angle(c);
        normal += phi*n;
    }
    return unit(normal); // placeholder
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    // TODO
    Vector3 normal{0, 0, 0};
    for (Corner c : v.adjacentCorners()) {
        Vector3 n = faceNormal(c.face());
        Vector3 eij = halfedgeVector(c.halfedge());
        Vector3 eik = halfedgeVector(c.halfedge().next().next().twin());
        // normal += norm(cross(eij, eik)/(norm2(eij)*norm2(eik))) * n;
        normal += cross(eij, eik)/(norm2(eij)*norm2(eik));
    }
    return unit(normal); // placeholder
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    // TODO
    Vector3 normal{0, 0, 0};
    for (Corner c : v.adjacentCorners()) {
        Vector3 n = faceNormal(c.face());
        double area = faceArea(c.face());
        normal += area * n;
    }
    return unit(normal); // placeholder
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    // TODO
    Vector3 normal{0, 0, 0};
    for (Halfedge he : v.outgoingHalfedges()) {
        double theta = dihedralAngle(he);
        normal += 0.5*theta*unit(halfedgeVector(he));
    }
    return -unit(normal); // placeholder
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    // TODO
    Vector3 normal{0, 0, 0};
    for (Halfedge he : v.outgoingHalfedges()) {
        normal += 0.5*(cotan(he)+cotan(he.twin())) * halfedgeVector(he);
    }
    return -unit(normal); // placeholder
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    // TODO
    double defect = 0.0;
    for (Corner c : v.adjacentCorners()) {
        defect += angle(c);
    }
    return 2*PI - defect; // placeholder
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    // TODO
    double total = 0.0;
    for (Vertex v : mesh.vertices()) {
        total += angleDefect(v);
    }
    return total; // placeholder
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    // TODO
    double curvature = 0.0;
    for (Halfedge he : v.outgoingHalfedges()) {
        curvature += dihedralAngle(he)*norm(halfedgeVector(he));
    }
    return 0.5*curvature; // placeholder
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    // TODO
    double dualArea = 0.0;
    for (Corner c : v.adjacentCorners()) {
        Halfedge eij = c.halfedge();
        Halfedge eik = c.halfedge().next().next().twin();
        dualArea += norm2(halfedgeVector(eij))*cotan(eij) + norm2(halfedgeVector(eik))*cotan(eik.twin());
    }
    return (1./8.)*dualArea; // placeholder
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
    double H = scalarMeanCurvature(v) / circumcentricDualArea(v);
    double K = angleDefect(v) / circumcentricDualArea(v);
    double sqrt = std::sqrt(H*H-K);
    return std::make_pair(H-sqrt, H+sqrt); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral