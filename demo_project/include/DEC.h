#ifndef _DEG_
#define _DEG_
#include "geometry.h"
#include"Eigen/Sparse"
#include"Eigen/Dense"

extern CMyMesh mesh;

using namespace MeshLib;
template <typename M, typename V, typename E, typename F, typename H>
class DEC {
public:

	static void testFun() {
		std::cout << "mesh " << mesh.vertices().size() << std::endl;
	}
	//DEC() {}
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static Eigen::SparseMatrix<double> buildHodgeStar0Form(Geometry<M, V, E, F, H> *geometry, std::map<int, int> vertexIndex) {
		//let vertices = geometry.mesh.vertices;
		int V_ = geometry->mesh->numVertices();
		std::vector<Eigen::Triplet<double>> tripletlist;
		Eigen::Triplet<double>* T = new Eigen::Triplet<double>(V_, V_);
		for (M::MeshVertexIterator mviter(geometry->mesh);!mviter.end();++mviter) {
			V* v = mviter.value();
			int i = vertexIndex[v->id()];
			double area = geometry->barycentricDualArea(v);

			tripletlist.push_back(Eigen::Triplet<double>(i, i, area));
		}

		Eigen::SparseMatrix<double> sp(V_, V_);
		sp.setFromTriplets(tripletlist.begin(), tripletlist.end());
		return sp;
	}

	///**
	// * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	// * @static
	// * @param {module:Core.Geometry} geometry The geometry of a mesh.
	// * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	// * @returns {module:LinearAlgebra.SparseMatrix}
	// */
	//static Eigen::SparseMatrix<double> buildHodgeStar1Form(Geometry<M, V, E, F, H>* geometry, std::map<int, int> edgeIndex) {
	//	int E_ = geometry->getMesh()->numEdges();
	//	std::cout << geometry->getMesh()->numVertices() << std::endl;
	//	std::vector<Eigen::Triplet<double>> tripletlist;
	//	for (M::MeshEdgeIterator meiter(geometry->getMesh()); !meiter.end(); ++meiter) {
	//		E* e = meiter.value();
	//		int i = edgeIndex[e->id()];
	//		double w = (geometry->cotan(static_cast<H*>(e->halfedge(0))) + geometry->cotan(static_cast<H*> (e->halfedge(1)))) / 2;

	//		tripletlist.push_back(Eigen::Triplet<double>(i, i, w));
	//	}

	//	Eigen::SparseMatrix<double> sp(E_, E_);
	//	sp.setFromTriplets(tripletlist.begin(), tripletlist.end());
	//	return sp;
	//}

	///**
	// * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	// * By convention, the area of a vertex is 1.
	// * @static
	// * @param {module:Core.Geometry} geometry The geometry of a mesh.
	// * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	// * @returns {module:LinearAlgebra.SparseMatrix}
	// */
	//static Eigen::SparseMatrix<double> buildHodgeStar2Form(Geometry<M, V, E, F, H>* geometry, std::map<int, int> faceIndex) {
	//	int F_ = geometry->getMesh()->numFaces();
	//	std::vector<Eigen::Triplet<double>> tripletlist;
	//	for (M::MeshFaceIterator mfiter(geometry->getMesh()); !mfiter.end(); ++mfiter) {
	//		F* f = mfiter.value();
	//		int i = faceIndex[f->id()];
	//		double area = geometry->area(f);

	//		tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0 / area));
	//	}

	//	Eigen::SparseMatrix<double> sp(F_, F_);
	//	sp.setFromTriplets(tripletlist.begin(), tripletlist.end());
	//	return sp;
	//}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static Eigen::SparseMatrix<double> buildExteriorDerivative0Form(Geometry<M, V, E, F, H>* geometry, std::map<int, int> edgeIndex, std::map<int, int> vertexIndex) {
		int E_ = geometry->mesh->numEdges();
		int V_ = geometry->mesh->numVertices();
		//let T = new Triplet(E, V);
		std::vector<Eigen::Triplet<double>> tripletlist;
		for (M::MeshEdgeIterator meiter(geometry->mesh); !meiter.end(); ++meiter) {
			E* e = meiter.value();
			int i = edgeIndex[e->id()];
			int j = vertexIndex[e->halfedge(0)->target()->id()];
			int k;
			if (e->halfedge(1))
				k = vertexIndex[e->halfedge(1)->target()->id()];
			else
				k = j;

			tripletlist.push_back(Eigen::Triplet<double>(i, j, 1));
			tripletlist.push_back(Eigen::Triplet<double>(i, k, -1));
		}

		Eigen::SparseMatrix<double> sp(E_, V_);
		sp.setFromTriplets(tripletlist.begin(), tripletlist.end());
		return sp;
	}

	///**
	// * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	// * @static
	// * @param {module:Core.Geometry} geometry The geometry of a mesh.
	// * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	// * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	// * @returns {module:LinearAlgebra.SparseMatrix}
	// */
	//static Eigen::SparseMatrix<double> buildExteriorDerivative1Form(Geometry<M, V, E, F, H>* geometry, std::map<int, int> faceIndex, std::map<int, int> edgeIndex) {
	//	int E_ = geometry->getMesh()->numEdges();
	//	int F_ = geometry->getMesh()->numFaces();
	//	//let T = new Triplet(F, E);
	//	std::vector<Eigen::Triplet<double>> tripletlist;
	//	for (M::MeshFaceIterator mfiter(geometry->getMesh()); !mfiter.end(); ++mfiter) {
	//		F* f = mfiter.value();
	//		int i = faceIndex[f->id()];

	//		for (M::FaceHalfedgeIterator fhiter(f); !fhiter.end(); ++fhiter) {
	//			H* h = fhiter.value();
	//			int j = edgeIndex[h->edge()->id()];
	//			double sign = h->edge()->halfedge(0) == h ? 1 : -1;

	//			tripletlist.push_back(Eigen::Triplet<double>(i, j, sign));
	//		}
	//	}

	//	Eigen::SparseMatrix<double> sp(F_, E_);
	//	sp.setFromTriplets(tripletlist.begin(), tripletlist.end());
	//	return sp;
	//}
};

#endif