#ifndef _MY_MESH_
#define _MY_MESH_


#include "Mesh\Vertex.h"
#include "Mesh\Edge.h"
#include "Mesh\Face.h"
#include "Mesh\HalfEdge.h"
#include "Mesh\BaseMesh.h"

#include "Mesh\boundary.h"
#include "Mesh\iterators.h"
#include "Parser\parser.h"

#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

namespace MeshLib
{
    class CMyVertex;
    class CMyEdge;
    class CMyFace;
    class CMyHalfEdge;

    /*! \brief CMyVertex class
     *
     *   Vertex class for this demo
     *   Trait : Vertex color
     */
    class CMyVertex : public CVertex
    {
    public:

        CPoint direction;
        double length;
        /*! constructor */
        CMyVertex() : m_rgb(229.0 / 255.0, 162.0 / 255.0, 141.0 / 255.0) { direction = CPoint(0, 0, 0); length = 0; };
        
        /*! read vertex attributes */
        void _from_string();

        /*! write vertex attributes */
        void _to_string();

        /*! vertex color */
        CPoint& rgb() { return m_rgb; };
    protected:
        /*! vertex color */
        CPoint m_rgb;
    };

    inline void CMyVertex::_from_string()
    {
        CParser parser(m_string);
        for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
        {
            CToken* token = *iter;
            if (token->m_key == "uv") //CPoint2
            {
                token->m_value >> m_uv;
            }
            if (token->m_key == "rgb") // CPoint
            {
                token->m_value >> m_rgb;
            }
        }
    }

    inline void CMyVertex::_to_string()
    {
        CParser parser(m_string);
        parser._removeToken("uv");

        parser._toString(m_string);
        std::stringstream iss;

        iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

        if (m_string.length() > 0)
        {
            m_string += " ";
        }
        m_string += iss.str();
    }
    
    /*! \brief CMyEdge class
     *
     *   Edge class for this demo
     *   Trait : Edge sharp
     */
    class CMyEdge : public CEdge
    {
    public:
        /*! constructor */
        CMyEdge() :m_sharp(false) {};

        /*! read edge attributes */
        void _from_string();

        /*! write edge attributes */
        void _to_string();

        /*! sharp edge */
        bool& sharp() { return m_sharp; };

        double length() {
            return (this->halfedge(0)->target()->point() -
                this->halfedge(0)->source()->point()).norm();
        }
    protected:
        /*! sharp edge */
        bool m_sharp;
    };

    inline void CMyEdge::_from_string()
    {
        CParser parser(m_string);
        for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
        {
            CToken* token = *iter;
            if (token->m_key == "sharp") // bool
            {
                m_sharp = true;
            }
        }
    }

    inline void CMyEdge::_to_string()
    {
        CParser parser(m_string);
        parser._removeToken("sharp");

        parser._toString(m_string);
        std::stringstream iss;

        if (m_sharp)
            iss << "sharp";

        if (m_string.length() > 0)
        {
            m_string += " ";
        }
        m_string += iss.str();
    }

    /*! \brief CMyFace class
     *
     *   Face class for this demo
     *   Trait : Face normal
     */
    class CMyFace : public CFace
    {
    public:
        /*! face normal */
        CPoint& normal() { return m_normal; };
    protected:
        /*! face normal */
        CPoint m_normal;
    };

    /*! \brief CMyHalfEdge class
     *
     *   HalfEdge class for this demo
     */
    class CMyHalfEdge : public CHalfEdge
    {
    };

    /*! \brief CMyMesh class
     *
     *	Mesh class for this demo
     *
     */
    template<typename V, typename E, typename F, typename H>
    class MyMesh : public CBaseMesh<V, E, F, H>
    {
    public:
        typedef V V;
        typedef E E;
        typedef F F;
        typedef H H;

        typedef CBoundary<V, E, F, H>					CBoundary;
        typedef CLoop<V, E, F, H>						CLoop;

        typedef MeshVertexIterator<V, E, F, H>			MeshVertexIterator;
        typedef MeshEdgeIterator<V, E, F, H>			MeshEdgeIterator;
        typedef MeshFaceIterator<V, E, F, H>			MeshFaceIterator;
        typedef MeshHalfEdgeIterator<V, E, F, H>		MeshHalfEdgeIterator;

        typedef VertexVertexIterator<V, E, F, H>		VertexVertexIterator;
        typedef VertexEdgeIterator<V, E, F, H>			VertexEdgeIterator;
        typedef VertexFaceIterator<V, E, F, H>			VertexFaceIterator;
        typedef VertexInHalfedgeIterator<V, E, F, H>	VertexInHalfedgeIterator;
        typedef VertexOutHalfedgeIterator<V, E, F, H>	VertexOutHalfedgeIterator;

        typedef FaceVertexIterator<V, E, F, H>			FaceVertexIterator;
        typedef FaceEdgeIterator<V, E, F, H>			FaceEdgeIterator;
        typedef FaceHalfedgeIterator<V, E, F, H>		FaceHalfedgeIterator;

        /*!
         *  \brief output the mesh information
         */
        void outputMeshInfo();
        
        /*!
         *  \brief a demo to show how to use iterators of MeshLib
         */
        void testIterator();
    };

    typedef MyMesh<CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> CMyMesh;

    template<typename V, typename E, typename F, typename H>
    void MeshLib::MyMesh<V, E, F, H>::outputMeshInfo()
    {
        int nv = this->numVertices();
        int ne = this->numEdges();
        int nf = this->numFaces();

        std::cout << "#V=" << nv << "  ";
        std::cout << "#E=" << ne << "  ";
        std::cout << "#F=" << nf << "  ";

        int euler_char = nv - ne + nf;
        std::cout << "Euler's characteristic=" << euler_char << "  ";

        CBoundary boundary(this);
        std::vector<CLoop*>& loops = boundary.loops();
        int nb = loops.size();

        int genus = (2 - (euler_char + nb)) / 2;
        std::cout << "genus=" << genus << std::endl;
    }

    template<typename V, typename E, typename F, typename H>
    void MyMesh<V, E, F, H>::testIterator()
    {
        for (MeshVertexIterator viter(this); !viter.end(); ++viter)
        {
            V* pV = *viter;
            // you can do something to the vertex here
            // ...

            for (VertexVertexIterator vviter(pV); !vviter.end(); ++vviter)
            {
                V* pW = *vviter;
                // you can do something to the neighboring vertices with CCW
                // ...
            }

            for (VertexEdgeIterator veiter(pV); !veiter.end(); ++veiter)
            {
                E* pE = *veiter;
                // you can do something to the neighboring edges with CCW
                // ...
            }

            for (VertexFaceIterator vfiter(pV); !vfiter.end(); ++vfiter)
            {
                F* pF = *vfiter;
                // you can do something to the neighboring faces with CCW
                // ...
            }

            for (VertexInHalfedgeIterator vhiter(this, pV); !vhiter.end(); ++vhiter)
            {
                H* pH = *vhiter;
                // you can do something to the incoming halfedges with CCW
                // ...
            }
        }

        for (MeshEdgeIterator eiter(this); !eiter.end(); ++eiter)
        {
            E* pE = *eiter;
            // you can do something to the edge here
            // ...
        }

        for (MeshFaceIterator fiter(this); !fiter.end(); ++fiter)
        {
            F* pF = *fiter;
            // you can do something to the face here
            // ...
        }

        //there are some other iterators which you can find them in class MyMesh

        std::cout << "Iterators test OK.\n";
    }
}

#endif // !_MY_MESH_
