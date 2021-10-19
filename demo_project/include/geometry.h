#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>
#include "MyMesh.h"

#define PI 3.14159


template<typename M, typename V, typename E, typename F, typename H>
class Geometry
{
public:
	M* mesh;

	double cotan(H* h);
	MeshLib::CPoint faceNormal(H* he);
	MeshLib::CPoint faceNormal(F* he);
	double dihedralAngle(H* h);

	std::map<int, int> getedgeIndex();
	std::map<int, int> getvertexIndex();
	double meanEdgeLength();

	double area(F* face) {
		CPoint p1 = face->halfedge()->target()->point() - face->halfedge()->source()->point();
		CPoint p2 = face->halfedge()->target()->point() - face->halfedge()->he_next()->target()->point();
		return (p1 ^ p2 / 2).norm();
	}

	/**
	 * Computes the barycentric dual area of a vertex.
	 * @method module:Core.Geometry#barycentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose barycentric dual area needs to be computed.
	 * @returns {number}
	 */
	double barycentricDualArea(V* v) {
		double area = 0.0;
		H* he_begin = (H*)v->most_ccw_out_halfedge();
		H* he_prev = he_begin;
		H* he = (H*)he_prev->he_sym()->he_next();
		do
		{
			MeshLib::CPoint p1 = he_prev->target()->point() - he_prev->source()->point();
			p1 /= p1.norm();
			MeshLib::CPoint p2 = he->target()->point() - he->source()->point();
			p2 /= p2.norm();

			area += std::acos(std::max(-1.0, std::min(1.0, p1 * p2)));

			he_prev = he;
			he = (H*)he_prev->he_sym()->he_next();

		} while (he_prev != he_begin);

		return area;
	}

	/**
	 * Computes the (integrated) scalar mean curvature at a vertex.
	 * @method module:Core.Geometry#scalarMeanCurvature
	 * @param {module:Core.Vertex} v The vertex whose mean curvature needs to be computed.
	 * @returns {number}
	 */
	double circumcentricDualArea(V* vertex) {
		double area = 0.0;
		H* he_begin = (H*)vertex->most_ccw_out_halfedge();
		H* he = he_begin;
		do
		{
			double a = (he->target()->point() - he->source()->point())
				^ (he->he_next()->target()->point() - he->he_next()->source()->point()) / 6;

			area += a;
			he = (H*)he->he_sym()->he_next();
		} while (he != he_begin);

		return area;
	}

	

	/**
	 * Computes the (integrated) scalar mean curvature at a vertex.
	 * @method module:Core.Geometry#scalarMeanCurvature
	 * @param {module:Core.Vertex} v The vertex whose mean curvature needs to be computed.
	 * @returns {number}
	 */
	double scalarMeanCurvature(V* vertex) {
		double sum = 0.0;
		H* he_begin = (H*)vertex->most_ccw_out_halfedge();
		H* he = he_begin;
		do
		{
			sum += 0.5 * ((E*)he->edge())->length() * dihedralAngle(he);
			he = (H*)he->he_sym()->he_next();
		} while (he!=he_begin);

		return sum;
	}

	double angleDefect(V* vertex) {
		double angleSum = 0.0;

		H* he_begin = (H*)vertex->most_ccw_out_halfedge();
		H* he_prev = he_begin;
		H* he = (H*)he_prev->he_sym()->he_next();
		do
		{
			MeshLib::CPoint p1 = he_prev->target()->point() - he_prev->source()->point();
			p1 /= p1.norm();
			MeshLib::CPoint p2 = he->target()->point() - he->source()->point();
			p2 /= p2.norm();

			angleSum += std::acos(std::max(-1.0, std::min(1.0, p1 * p2)));

			he_prev = he;
			he = (H*)he_prev->he_sym()->he_next();

		} while (he_prev != he_begin);


		return(vertex->boundary() ? PI - angleSum : 2 * PI - angleSum);
	}

	double scalarGaussCurvature(V* vertex) {
		return angleDefect(vertex);
	}
	 

};

template<typename M, typename V, typename E, typename F, typename H>
double Geometry<M, V, E, F, H>::cotan(H* h)
{
	if (h->edge()->boundary()) {
		return 0.0;
	}
	MeshLib::CPoint v1 = h->source()->point();
	MeshLib::CPoint v0 = h->he_prev()->source()->point();
	MeshLib::CPoint v2 = h->he_next()->source()->point();
	MeshLib::CPoint u = v1 - v0;
	MeshLib::CPoint v = v2 - v0;

	double udotv = u * v;
	MeshLib::CPoint ucrossv = u ^ v;
	double ans = udotv / ucrossv.norm();
	return ans;
}

template<typename M, typename V, typename E, typename F, typename H>
MeshLib::CPoint Geometry<M, V, E, F, H>::faceNormal(H* he)
{
	MeshLib::CPoint p = (he->target()->point() - he->source()->point()) ^ (he->he_next()->source()->point() - he->he_next()->target()->point());
	p /= p.norm();
	return(p);
	
}

template<typename M, typename V, typename E, typename F, typename H>
MeshLib::CPoint Geometry<M, V, E, F, H>::faceNormal(F* face)
{
	return(faceNormal((H*)face->halfedge()));

}

//template<typename M, typename V, typename E, typename F, typename H>
//double Geometry<M, V, E, F, H>::dihedralAngle(H* he) {
//	if (he->edge()->boundary()) return 0.0;
//	H* he_an;
//	MeshLib::CPoint faceNormal1 = (he->target()->point() - he->source()->point()) ^ (he->he_next()->source()->point() - he->he_next()->target()->point());
//	faceNormal1 /= faceNormal1.norm();
//	he_an = (H*)he->he_sym();
//	MeshLib::CPoint faceNormal2 = (he_an->target()->point() - he_an->source()->point()) ^ (he_an->he_next()->source()->point() - he_an->he_next()->target()->point());
//	faceNormal2 /= faceNormal2.norm();
//	MeshLib::CPoint w = he->target()->point() - he->source()->point();
//	w /= w.norm();
//
//	double cosTheta = faceNormal1 * faceNormal2;
//	double sinTheta = (faceNormal1 ^ faceNormal2) * w;
//
//	return atan2(sinTheta, cosTheta);
//}

#include <Eigen/Dense>
#include <Eigen/Sparse>

template<typename M, typename V, typename E, typename F, typename H>
double Geometry<M, V, E, F, H>::dihedralAngle(H* he)
{
	if (he->edge()->boundary()) {
		return 0.0;
	}
	MeshLib::CPoint n1 = faceNormal(he);
	MeshLib::CPoint n2 = faceNormal((H*)(he->he_sym()));
	MeshLib::CPoint w = he->target()->point() - he->source()->point();
	w /= w.norm();

	double cosTheta = n1 * n2;
	double sinTheta = (n1 ^ n2) * w;

	return atan2(sinTheta, cosTheta);
}

template<typename M, typename V, typename E, typename F, typename H>
inline std::map<int, int> Geometry<M, V, E, F, H>::getedgeIndex()
{
	int i = 0;
	std::map<int, int> edgeIndex;
	for (M::MeshEdgeIterator meiter(this->mesh); !meiter.end(); ++meiter) {
		E* e = meiter.value();
		edgeIndex[e->id()] = i;
		i++;
	}
	return edgeIndex;
}

template<typename M, typename V, typename E, typename F, typename H>
inline std::map<int, int> Geometry<M, V, E, F, H>::getvertexIndex()
{
	int i = 0;
	std::map<int, int> vertexIndex;
	for (M::MeshVertexIterator mviter(this->mesh); !mviter.end(); ++mviter) {
		V* v = mviter.value();
		vertexIndex[v->id()] = i;
		i++;
	}
	return vertexIndex;
}

template<typename M, typename V, typename E, typename F, typename H>
inline double Geometry<M, V, E, F, H>::meanEdgeLength()
{
	if (this->mesh == NULL)
		return -1;


	list<E*>::iterator iter = mesh->edges().begin();
	double sumLength = 0;
	while (iter != mesh->edges().end())
	{
		sumLength += (*iter)->length();
		iter++;
	}
	sumLength /= mesh->edges().size();
	return sumLength;
}

#endif // !GEOMETRY_H