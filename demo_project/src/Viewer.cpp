#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include "viewer/Arcball.h"                           /*  Arc Ball  Interface         */
#include "bmp/RgbImage.h"
#include "MyMesh.h"
#include "geometry.h"
#include "colormap.h"
#include "DEC.h"
#include "arrowPoint.h"

#include <Eigen/Dense>q

using namespace MeshLib;
using namespace Eigen;
using namespace std;

/* window width and height */
int win_width, win_height;
int gButton;
int startx, starty;
int shadeFlag = 1;
bool showMesh = true;
bool showArrows = false;
bool showUV = false;
bool showAxis = false;

/* rotation quaternion and translation vector for the object */
CQrot       ObjRot(0, 0, 1, 0);
CPoint      ObjTrans(0, 0, 0);

/* global mesh */
CMyMesh mesh;

map<int, ArrowPoint> arrowMesh;

/* arcball object */
CArcball arcball;

int textureFlag = 2;
/* texture id and image */
GLuint texName;
RgbImage image;
bool hasTexture = false;

enum class FormName
{
    Primal_0,
    Primal_1,
    Primal_2,
    Dual_0,
    Dual_1,
    Dual_2
};
FormName currentFormName = FormName::Primal_0;
Matrix<double, Dynamic, Dynamic> currentForm;
Geometry<CMyMesh, CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> geometry;

double gaussian(double x, double a, double b) {
    return a * exp(-(x * x) / (b * b));
}

void initArrowMesh() {
    for (list<CMyFace*>::iterator iter = mesh.faces().begin(); iter != mesh.faces().end(); iter++) {
        CMyFace* face = *iter;

        CMyHalfEdge* he = (CMyHalfEdge*)face->halfedge();
        CPoint p = he->target()->point() + he->source()->point() + he->he_next()->target()->point();
        p /= 3;

        ArrowPoint arrowPoint(p);
        arrowMesh.insert(make_pair(face->id(), arrowPoint));
    }
}

Matrix<double, Dynamic, Dynamic> generateRandomFormOnVertices(bool scaleByArea = false) {
    int v = mesh.vertices().size();
    vector<CPoint*> c;

    srand(time(0));
    double perks = max(2, rand() % 10);
    
    for (int i = 0; i < perks; i++) {
        int t = rand() % v;
        list<CMyVertex*> vertices = mesh.vertices();
        list<CMyVertex*>::iterator iter = vertices.begin();
        for (int j = 0; j < t; j++) {
            iter++;
        }
        CMyVertex* vertex = *iter;
        c.push_back(&vertex->point());
    }
    Matrix<double, Dynamic, Dynamic> form;
    form.resize(v, 1);

    list<CMyVertex*> vertices = mesh.vertices();
    for (list<CMyVertex*>::iterator iter = vertices.begin();
        iter != vertices.end(); iter++) {
        CMyVertex* vertex = *iter;
        int i = vertex->id();
        CPoint p = vertex->point();

        double sum = 0;
        for (int j = 0; j < perks; j++) {
            sum += gaussian((*(c[j]) - p).norm(), 1.0, 0.5);
        }
        if (scaleByArea) sum *= geometry.barycentricDualArea(vertex);
        form(i - 1, 0) = sum;
    }
    return form;
}

void updateColorsPrimalForm(Matrix<double, Dynamic, Dynamic>& form, bool is0Form) {

    double maxForm = form.maxCoeff();
    double minForm = form.minCoeff();

    if (currentFormName == FormName::Dual_1 || currentFormName == FormName::Primal_1) {
        if (maxForm <= 1e-3) mapToColor(0, minForm, maxForm, "coolwarm");
        else {
            list<CMyVertex*> vertices = mesh.vertices();
            for (list<CMyVertex*>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
                (*iter)->rgb()[0]= (*iter)->rgb()[1]= (*iter)->rgb()[2]=1;
            }
        }
    } else {
        for (list<CMyFace*>::iterator iter = mesh.faces().begin(); iter != mesh.faces().end(); iter++) {
            CMyFace* face = *iter;
            double A = geometry.area(*iter);
            CPoint color;
            if (!is0Form) {
                color = mapToColor(form(face->id(), 0) / A , minForm, maxForm, "coolwarm");
            }

            CMyHalfEdge* he = (CMyHalfEdge*)face->halfedge();
            do
            {
                CMyVertex* vertice = (CMyVertex*)he->target();

                if (is0Form) vertice->rgb() = mapToColor(form(vertice->id()-1, 0), minForm, maxForm, "coolwarm");
                else vertice->rgb() = color;
                he = (CMyHalfEdge*)he->he_next();
            } while (he!= (CMyHalfEdge*)face->halfedge());
        }
    }
}

void updateColorsDualForm(Matrix<double, Dynamic, Dynamic>& form, bool is2Form) {

    double maxForm = form.maxCoeff();
    double minForm = form.minCoeff();

    if (currentFormName == FormName::Dual_1 || currentFormName == FormName::Primal_1) {
        if (maxForm <= 1e-3) mapToColor(0, minForm, maxForm, "coolwarm");
        else {
            list<CMyVertex*> vertices = mesh.vertices();
            for (list<CMyVertex*>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
                (*iter)->rgb()[0] = (*iter)->rgb()[1] = (*iter)->rgb()[2] = 1;
            }
        }
    }
    else {
        int k = 0;
        for (list<CMyFace*>::iterator iter = mesh.faces().begin(); iter != mesh.faces().end(); iter++) {
            CMyFace* f = *iter;

            CMyHalfEdge* he = (CMyHalfEdge*)f->halfedge();
            for (int i = 0; i < 3; i++) {

                CMyVertex* v = (CMyVertex*)he->target();
                CPoint color;
                double avg1RingFormValue = 0;

                if (is2Form) {
                    double A = geometry.barycentricDualArea(v);
                    color = mapToColor(form(v->id(), 0) / A, minForm, maxForm, "coolwarm");
                }
                else {
                    int n = 0;
                    CMyHalfEdge* he_begin = (CMyHalfEdge*)v->most_ccw_out_halfedge();
                    CMyHalfEdge* he = (CMyHalfEdge*)he_begin->he_sym()->he_next();
                    do
                    {
                        avg1RingFormValue += form(f->id() - 1, 0);
                        he = (CMyHalfEdge*)he->he_sym()->he_next();
                        n++;

                    } while (he != he_begin);

                    avg1RingFormValue /= n;
                }

                CMyHalfEdge* h = (CMyHalfEdge*)he->he_prev();
                if (h->edge()->boundary()) {
                    for (int i = 0; i < 6; i++) {
                        if (!is2Form) {
                            double formValue = form(f->id() - 1, 0);
                            if (i == 0 || i == 3)formValue = avg1RingFormValue;
                            else if (i == 4)formValue = form(h->he_prev()->he_sym()->face()->id() - 1, 0);
                            else if (i == 5)formValue = form(h->he_prev()->he_sym()->face()->id() - 1, 0);

                            color = mapToColor(formValue, minForm, maxForm, "coolwarm");
                        }
                    }

                    
                }




                he = (CMyHalfEdge*)he->he_next();
            }

        }

    }
}

void drawArrow( CPoint& point, const CPoint& direction) {
    double length = 2;
    glBegin(GL_LINES);
    glColor3f(0.0, 0.0, 0.0);	//red
    //CPoint p = direction / direction.norm();
    glVertex3d(point[0] - 0.5 * direction[0] * length, point[1] - 0.5 * direction[1] * length, point[2] - 0.5 * direction[2] * length);
    glVertex3d(point[0] + 0.5 * direction[0] * length, point[1] + 0.5 * direction[1] * length, point[2] + 0.5 * direction[2] * length);
    glEnd();

    glBegin(GL_TRIANGLES); 
    double arrowSize = 0.1;
    glVertex3f(point[0] + (0.5 * direction[0]) * (length + arrowSize) , point[1] + (0.5 * direction[1]) * (length + arrowSize), point[2] + (0.5 * direction[2]) * (length + arrowSize));
    glVertex3d(point[0] + 0.5 * direction[0] * length - 2 * arrowSize * direction[1] * length, point[1] + 0.5 * direction[1] * length - 2 * arrowSize * direction[0] * length, point[2] + 0.5 * direction[2] * length);
    glVertex3d(point[0] + 0.5 * direction[0] * length + 2 * arrowSize * direction[1] * length, point[1] + 0.5 * direction[1] * length + 2 * arrowSize * direction[0] * length, point[2] + 0.5 * direction[2] * length);
    glEnd();
}

void drawArrows() {
    for (map<int, ArrowPoint>::iterator iter = arrowMesh.begin(); iter != arrowMesh.end(); iter++) {
        pair<int, ArrowPoint> p = *iter;

        drawArrow(p.second.point, p.second.direction);
    }
    
    /*for (list<CMyVertex*>::iterator iter = mesh.vertices().begin(); iter != mesh.vertices().end(); iter++) {
        drawArrow((*iter)->point(), CPoint(1, 2, 2));
    }*/
}

map<int, CPoint> interpolateWhitney(Matrix<double, Dynamic, Dynamic>& form) {
    map<int, CPoint> field;
    for (list<CMyFace*>::iterator iter = mesh.faces().begin(); iter != mesh.faces().end(); iter++) {
        CMyFace* face = *iter;
        CMyHalfEdge* he = (CMyHalfEdge*)face->halfedge();

        CPoint pi = he->target()->point();
        CPoint pj = he->he_next()->target()->point();
        CPoint pk = he->source()->point();
        CPoint eij = pj - pi;
        CPoint ejk = pk - pj;
        CPoint eki = pi - pk;
        //cout << "id  " << he->edge()->id() << endl;
        double cij = form(he->edge()->id() -1 , 0);
        double cjk = form(he->he_next()->edge()->id() -1 , 0);
        double cki = form(he->he_prev()->edge()->id() -1 , 0);

        if (he->edge()->halfedge(0) != he)cij *= -1;
        if (he->he_next()->edge()->halfedge(0) != he->he_next())cjk *= -1;
        if (he->he_prev()->edge()->halfedge(0) != he->he_prev())cki *= -1;

        double A = geometry.area(face);
        CPoint N = face->normal();
        
        CPoint a = (eki - ejk) * cij;
        CPoint b = (eij - eki) * cjk;
        CPoint c = (ejk - eij) * cki;

        CPoint res = N ^ ((a + b + c) / (A * 6));
        field.insert(make_pair(face->id(), res));
    }
    return field;
}

void updatePrimal1FormMesh() {
    cout << "updatePrimal1FormMesh" << endl;
    // interpolate 1 form to a face field
    map<int, CPoint> primal1FormField = interpolateWhitney(currentForm);
    cout << "size  " << primal1FormField.size() << endl;
    double length = 0.3 * geometry.meanEdgeLength();

    for (list<CMyFace*>::iterator iter = mesh.faces().begin(); iter != mesh.faces().end(); iter++) {
        CMyFace* face = *iter;
        

        CPoint N = geometry.faceNormal(face);
        map<int, CPoint>::iterator p = primal1FormField.find(face->id());
        CPoint field = (*p).second * length;

        double norm = field.norm();
        if (norm > length) {
            field *= length / norm;
        }
        
        (*arrowMesh.find(face->id())).second.direction = field;

       /* cout << face->id() << " " << (*arrowMesh.find(face->id())).second.direction(0)
            << " " << (*arrowMesh.find(face->id())).second.direction(1)
            << " " << (*arrowMesh.find(face->id())).second.direction(2) << endl;*/
    }

    // update positions
    /*let positions = primal1FormMesh.geometry.attributes.position.array;
    for (let f of mesh.faces) {
        let C = geometry.circumcenter(f);
        let N = geometry.faceNormal(f);

        let field = primal1FormField[f].times(length);
        clampFieldLength(field, length);
        setArrow(positions, f.index, C.minus(field), C.plus(field), N);
    }

    primal1FormMesh.geometry.attributes.position.needsUpdate = true;*/
}




void updateFormViz() {
    if (currentFormName == FormName::Primal_0 || currentFormName == FormName::Primal_1
        || currentFormName == FormName::Primal_2) {

        if (currentFormName == FormName::Primal_1) {
            updateColorsPrimalForm(currentForm, true);
            updatePrimal1FormMesh();
            showArrows = true;
        } else {
            showArrows = false;
            if (currentFormName == FormName::Primal_0) {
                updateColorsPrimalForm(currentForm, true);
            }
            else {
                updateColorsPrimalForm(currentForm, false);
            }
        }
    }
    else {
        if (currentFormName == FormName::Dual_1) {
            updateColorsDualForm(currentForm, false);
            updatePrimal1FormMesh();
            showArrows = true;
        }
        else {
            showArrows = false;
            if (currentFormName == FormName::Dual_0) {
                updateColorsDualForm(currentForm, false);
            }
            else {
                updateColorsDualForm(currentForm, true);
            }
        }
    }
}

void functionD() {

    Eigen::SparseMatrix<double> d;
    cout << "begin d" << endl;

    switch (currentFormName)
    {
    case FormName::Primal_0:
        d = DEC<CMyMesh, CMyVertex, CMyEdge, CMyFace, CMyHalfEdge>::
            buildExteriorDerivative0Form(&geometry, geometry.getedgeIndex(), geometry.getvertexIndex());
        currentFormName = FormName::Primal_1;
        break;
    case FormName::Primal_1:
        currentFormName = FormName::Primal_2;
        break;
    case FormName::Dual_0:
        currentFormName = FormName::Dual_1;
        break;
    case FormName::Dual_1:
        currentFormName = FormName::Dual_2;
        break;
    default:
        break;
    }

    currentForm = d * currentForm;
    cout << "currentForm " << currentForm.cols() << " " << currentForm.rows() << endl;
    updateFormViz();
    cout << "end d" << endl;
}

void functionStar() {

    Eigen::SparseMatrix<double> star;
    cout << "begin star" << endl;

    switch (currentFormName)
    {
    case FormName::Primal_0:
        star = DEC<CMyMesh, CMyVertex, CMyEdge, CMyFace, CMyHalfEdge>::
            buildHodgeStar0Form(&geometry, geometry.getvertexIndex());
        currentFormName = FormName::Dual_2;
        break;
    case FormName::Primal_1:
        currentFormName = FormName::Dual_1;
        break;
    case FormName::Primal_2:
        currentFormName = FormName::Dual_0;
        break;
    default:
        break;
    }


    currentForm = star * currentForm;
    cout << "currentForm " << currentForm.cols() << " " << currentForm.rows() << endl;
    updateFormViz();
    cout << "end star" << endl;
}


void randomize() {
    if (currentFormName == FormName::Primal_0 || currentFormName == FormName::Dual_2) {
        currentForm = generateRandomFormOnVertices(currentFormName == FormName::Dual_2);
    }

    updateFormViz();
}

//copy frame buffer to an image
/*! save frame buffer to an image "snap_k.bmp"
*/
void readFrameBuffer()
{
    static int id = 0;
    GLfloat* buffer = new GLfloat[win_width * win_height * 3];
    assert(buffer);
    glReadBuffer(GL_FRONT_LEFT);
    glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_FLOAT, buffer);

    RgbImage image(win_height, win_width);

    for (int i = 0; i < win_height; i++)
        for (int j = 0; j < win_width; j++)
        {
            float r = buffer[(i * win_width + j) * 3 + 0];
            float g = buffer[(i * win_width + j) * 3 + 1];
            float b = buffer[(i * win_width + j) * 3 + 2];

            image.SetRgbPixelf(i, j, r, g, b);
        }
    delete[]buffer;

    char name[256];
    std::ostringstream os(name);
    os << "snape_" << id++ << ".bmp";
    image.WriteBmpFile(os.str().c_str());

}

/*! initialize bitmap image texture */
void initializeBmpTexture()
{
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,   GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    //	int ImageWidth  = image.GetNumRows();
    //	int ImageHeight = image.GetNumCols();
    int ImageWidth = image.GetNumCols();
    int ImageHeight = image.GetNumRows();
    GLubyte* ptr = (GLubyte*)image.ImageData();

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
        ImageWidth,
        ImageHeight,
        0,
        GL_RGB,
        GL_UNSIGNED_BYTE,
        ptr);

    if (textureFlag == 1)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    else if (textureFlag == 2)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_TEXTURE_2D);
}

/*! setup the object, transform from the world to the object coordinate system */
void setupObject(void)
{
    double rot[16];

    glTranslated(ObjTrans[0], ObjTrans[1], ObjTrans[2]);
    ObjRot.convert(rot);
    glMultMatrixd((GLdouble*)rot);
}

/*! the eye is always fixed at world z = +5 */
void setupEye(void) {
    glLoadIdentity();
    gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
}

/*! setup light */
void setupLight()
{
    GLfloat lightOnePosition[4] = { 0, 0, 1, 0 };
    GLfloat lightTwoPosition[4] = { 0, 0, -1, 0 };
    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
    glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);
}

void renderStrokeFontString(
    float x,
    float y,
    float z,
    void* font,
    char* string) {

    char* c;
    glPushMatrix();
    glTranslatef(x, y, z);
    glScalef(1 / 1800.0, 1 / 1800.0, 1 / 1800.0);

    for (c = string; *c != '\0'; c++) {
        glutStrokeCharacter(font, *c);
    }

    glPopMatrix();
}

/*! draw axis */
void drawAxis()
{
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);

    glLineWidth(1.5);
    //x axis
    glColor3f(1.0, 0.0, 0.0);	//red
    glBegin(GL_LINES);
    glVertex3d(0, 0, 0);
    glVertex3d(1, 0, 0);
    glEnd();
    renderStrokeFontString(1, 0, 0, GLUT_STROKE_MONO_ROMAN, "X");

    //y axis
    glColor3f(0.0, 1.0, 0);		//green
    glBegin(GL_LINES);
    glVertex3d(0, 0, 0);
    glVertex3d(0, 1, 0);
    glEnd();
    renderStrokeFontString(0, 1, 0, GLUT_STROKE_MONO_ROMAN, "Y");

    //z axis
    glColor3f(0.0, 0.0, 1.0);	//blue
    glBegin(GL_LINES);
    glVertex3d(0, 0, 0);
    glVertex3d(0, 0, 1);
    glEnd();
    renderStrokeFontString(0, 0, 1, GLUT_STROKE_MONO_ROMAN, "Z");

    glLineWidth(1.0);
}

/*! draw mesh */
void drawMesh()
{
    glEnable(GL_LIGHTING);
    if (hasTexture)
        glBindTexture(GL_TEXTURE_2D, texName);

    for (CMyMesh::MeshFaceIterator fiter(&mesh); !fiter.end(); ++fiter)
    {
        glBegin(GL_POLYGON);
        CMyFace* pf = *fiter;
        for (CMyMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
        {
            CMyVertex* v = *fviter;
            CPoint& pt = v->point();
            CPoint2& uv = v->uv();
            CPoint& rgb = v->rgb();
            CPoint n;
            switch (shadeFlag)
            {
            case 0:
                n = pf->normal();
                break;
            case 1:
                n = v->normal();
                break;
            }
            glNormal3d(n[0], n[1], n[2]);
            glTexCoord2d(uv[0], uv[1]);
            glColor3f(rgb[0], rgb[1], rgb[2]);
            glVertex3d(pt[0], pt[1], pt[2]);
        }
        glEnd();
    }
}

/*! draw uv */
void drawUv()
{
    glEnable(GL_LIGHTING);
    if (hasTexture)
        glBindTexture(GL_TEXTURE_2D, texName);

    for (CMyMesh::MeshFaceIterator fiter(&mesh); !fiter.end(); ++fiter)
    {
        glBegin(GL_POLYGON);
        CMyFace* pf = *fiter;
        for (CMyMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
        {
            CMyVertex* v = *fviter;
            CPoint2& uv = v->uv();
            CPoint& rgb = v->rgb();
            CPoint n;
            switch (shadeFlag)
            {
            case 0:
                n = pf->normal();
                break;
            case 1:
                n = v->normal();
                break;
            }
            glNormal3d(n[0], n[1], n[2]);
            glTexCoord2d(uv[0], uv[1]);
            glColor3f(rgb[0], rgb[1], rgb[2]);
            glVertex3d(uv[0], uv[1], 0);
        }
        glEnd();
    }
}

void drawSharpEdges()
{
    glDisable(GL_LIGHTING);

    glLineWidth(1.5);
    glColor3f(1, 0, 0);
    glBegin(GL_LINES);
    for (CMyMesh::MeshEdgeIterator eiter(&mesh); !eiter.end(); ++eiter)
    {
        CMyEdge* pE = *eiter;
        if (pE->sharp() == true)
        {
            CMyVertex* p0 = mesh.edgeVertex1(pE);
            CMyVertex* p1 = mesh.edgeVertex2(pE);
            glColor3f(1.0f, 0.0f, 0.0f);
            glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
            glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
        }
    }
    glEnd();
    glLineWidth(1.0);
}

/*! display call back function
*/
void display()
{
    /* clear frame buffer */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    setupLight();
    /* transform from the eye coordinate system to the world system */
    setupEye();
    glPushMatrix();
    /* transform from the world to the ojbect coordinate system */
    setupObject();

    /* draw the axis */
    if (showAxis)
        drawAxis();

    /* draw sharp edges */
    drawSharpEdges();

    if(showArrows)
        drawArrows();

    /* draw the mesh */
    switch (textureFlag)
    {
    case 0:
        glDisable(GL_TEXTURE_2D);
        break;
    case 1:
        glEnable(GL_TEXTURE_2D);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        break;
    case 2:
        glEnable(GL_TEXTURE_2D);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        break;
    }
    if (showMesh)
        drawMesh();
    if (showUV)
        drawUv();

    glPopMatrix();
    glutSwapBuffers();
}

/*! Called when a "resize" event is received by the window. */
void reshape(int w, int h)
{
    float ar;
    //std::cout << "w:" << w << "\th:" << h << std::endl;
    win_width = w;
    win_height = h;

    ar = (float)(w) / h;
    glViewport(0, 0, w, h);               /* Set Viewport */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // magic imageing commands
    gluPerspective(40.0, /* field of view in degrees */
        ar, /* aspect ratio */
        0.1, /* Z near */
        100.0 /* Z far */);

    glMatrixMode(GL_MODELVIEW);

    glutPostRedisplay();
}

/*! helper function to remind the user about commands, hot keys */
void help()
{
    printf("0  -  Show the coordinate axis\n");
    printf("1  -  Show the original mesh\n");
    printf("2  -  Show the uv parametrization\n");
    printf("w  -  Wireframe Display\n");
    printf("f  -  Flat Shading \n");
    printf("s  -  Smooth Shading\n");
    printf("?  -  Help Information\n");
    printf("esc - quit\n");
}

void specialKey(GLint key, GLint x, GLint y)
{
    glutPostRedisplay();
}

/*! Keyboard call back function */
void keyBoard(unsigned char key, int x, int y)
{
    switch (key)
    {
    case '0':
        showAxis = !showAxis;
        break;
    case '1':
        showMesh = !showMesh;
        break;
    case '2':
        showUV = !showUV;
        break;
    case '8':
        break;
    case 'd':
        functionD();
        break;
    case 'f':
        //Flat Shading
        glPolygonMode(GL_FRONT, GL_FILL);
        shadeFlag = 0;
        break;
    case 's':
        //Smooth Shading
        glPolygonMode(GL_FRONT, GL_FILL);
        shadeFlag = 1;
        break;
    case 'w':
        //Wireframe mode
        glPolygonMode(GL_FRONT, GL_LINE);
        break;
    case 't':
        textureFlag = (textureFlag + 1) % 3;
        break;
    case 'o':
        readFrameBuffer();
        break;
    case '?':
        help();
        break;
    case 27:
        exit(0);
        break;
    }
    glutPostRedisplay();
}

/*! setup GL states */
void setupGLstate() {
    GLfloat lightOneColor[] = { 1, 1, 1, 1.0 };
    GLfloat globalAmb[] = { .1, .1, .1, 1 };
    GLfloat lightOnePosition[] = { .0, 0.0, 1.0, 1.0 };
    GLfloat lightTwoPosition[] = { .0, 0.0, -1.0, 1.0 };

    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.35, 0.53, 0.70, 0);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);

    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, lightOneColor);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
    glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);

    const GLfloat specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 64.0f);

    GLfloat mat_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat mat_diffuse[] = { 0.01f, 0.01f, 0.01f, 1.0f };
    GLfloat mat_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat mat_shininess[] = { 32 };

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
}

/*! mouse click call back function */
void  mouseClick(int button, int state, int x, int y) {
    /* set up an arcball around the Eye's center
    switch y coordinates to right handed system  */

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        gButton = GLUT_LEFT_BUTTON;
        arcball = CArcball(win_width, win_height, x - win_width / 2, win_height - y - win_height / 2);
    }

    if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
        startx = x;
        starty = y;
        gButton = GLUT_MIDDLE_BUTTON;
    }

    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
        startx = x;
        starty = y;
        gButton = GLUT_RIGHT_BUTTON;
    }
    return;
}

/*! mouse motion call back function */
void mouseMove(int x, int y)
{
    CPoint trans;
    CQrot  rot;

    /* rotation, call arcball */
    if (gButton == GLUT_LEFT_BUTTON)
    {
        rot = arcball.update(x - win_width / 2, win_height - y - win_height / 2);
        ObjRot = rot * ObjRot;
        glutPostRedisplay();
    }

    /*xy translation */
    if (gButton == GLUT_MIDDLE_BUTTON)
    {
        double scale = 10. / win_height;
        trans = CPoint(scale * (x - startx), scale * (starty - y), 0);
        startx = x;
        starty = y;
        ObjTrans = ObjTrans + trans;
        glutPostRedisplay();
    }

    /* zoom in and out */
    if (gButton == GLUT_RIGHT_BUTTON) {
        double scale = 10. / win_height;
        trans = CPoint(0, 0, scale * (starty - y));
        startx = x;
        starty = y;
        ObjTrans = ObjTrans + trans;
        glutPostRedisplay();
    }

}


/*! Normalize mesh
* \param pMesh the input mesh
*/
void normalizeMesh(CMyMesh* pMesh)
{
    CPoint s(0, 0, 0);
    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        s = s + v->point();
    }
    s = s / pMesh->numVertices();

    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        CPoint p = v->point();
        p = p - s;
        v->point() = p;
    }

    double d = 0;
    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        CPoint p = v->point();
        for (int k = 0; k < 3; k++)
        {
            d = (d > fabs(p[k])) ? d : fabs(p[k]);
        }
    }

    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        CPoint p = v->point();
        p = p / d;
        v->point() = p;
    }
};

/*! Compute the face normal and vertex normal
* \param pMesh the input mesh
*/
void computeNormal(CMyMesh* pMesh)
{
    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        CPoint n(0, 0, 0);
        for (CMyMesh::VertexFaceIterator vfiter(v); !vfiter.end(); ++vfiter)
        {
            CMyFace* pF = *vfiter;

            CPoint p[3];
            CHalfEdge* he = pF->halfedge();
            for (int k = 0; k < 3; k++)
            {
                p[k] = he->target()->point();
                he = he->he_next();
            }

            CPoint fn = (p[1] - p[0]) ^ (p[2] - p[0]);
            pF->normal() = fn / fn.norm();
            n += fn;
        }

        n = n / n.norm();
        v->normal() = n;
    }
};

void initOpenGL(int argc, char* argv[])
{
    if (hasTexture)
        image.LoadBmpFile(argv[2]);

    /* glut stuff */
    glutInit(&argc, argv);                /* Initialize GLUT */
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
    glutCreateWindow("Mesh Viewer");	  /* Create window with given title */
    glViewport(0, 0, 600, 600);

    glutDisplayFunc(display);             /* Set-up callback functions */
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMove);
    glutKeyboardFunc(keyBoard);
    glutSpecialFunc(&specialKey);
    setupGLstate();

    if (hasTexture)
        initializeBmpTexture();

    glutMainLoop();                       /* Start GLUT event-processing loop */
}

/*! main function for viewer
*/
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: input.m [texture.bmp]" << std::endl;
        return -1;
    }

    if (argc > 2)
        hasTexture = true;

    std::string mesh_name(argv[1]);
    if (strutil::endsWith(mesh_name, ".obj"))
    {
        mesh.read_obj(mesh_name.c_str());
    }
    if (strutil::endsWith(mesh_name, ".m"))
    {
        mesh.read_m(mesh_name.c_str());
    }
    if (strutil::endsWith(mesh_name, ".off"))
    {
        mesh.read_off(mesh_name.c_str());
    }
    geometry.mesh = &mesh;
    cout<<"mean Length "<<geometry.meanEdgeLength()<<endl;

    normalizeMesh(&mesh);
    computeNormal(&mesh);
    initArrowMesh();

    //for (map<int, ArrowPoint>::iterator iter = arrowMesh.begin(); iter != arrowMesh.end(); iter++) {
    //    pair<int, ArrowPoint> pair = (*iter);
    //    cout << "i" << pair.first << "  (" << pair.second.point[0] << ", " << pair.second.point[1] << ", " << pair.second.point[2] << endl;
    //}

    randomize();

    mesh.outputMeshInfo();
    mesh.testIterator();

    initOpenGL(argc, argv);
    return 0;
}

