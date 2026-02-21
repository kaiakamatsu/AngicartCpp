#ifndef WX_H
#define WX_H 1

/* Converts (`x`, `y`, `z`) to a normalized coordinate (`nx`, `ny`, `nz`).
 * @nx the *x* coordinate of interest
 * @ny the *y* coordinate of interest
 * @nz the *z* coordinate of interest
 * @nx the normalized *x* coordinate
 * @ny the normalized *y* coordinate
 * @nz the normalized *z* coordinate
 *
 * Converts (`x`, `y`, `z`) to a normalized coordinate (`nx`, `ny`, `nz`) according to `voldimsGlobal`.
 */
template <class T>
void normalizedCoord(T x, T y, T z, double &nx, double &ny, double &nz){
	double largestExtent(voldimsGlobal[0]*voxdimsGlobal[0]);
	if(largestExtent < voldimsGlobal[1]*voxdimsGlobal[1])
		largestExtent = voldimsGlobal[1]*voxdimsGlobal[1];
	if(largestExtent < voldimsGlobal[2]*voxdimsGlobal[2])
		largestExtent = voldimsGlobal[2]*voxdimsGlobal[2];
	double extentFactors[] = {voldimsGlobal[0]*voxdimsGlobal[0]/largestExtent, voldimsGlobal[1]*voxdimsGlobal[1]/largestExtent, voldimsGlobal[2]*voxdimsGlobal[2]/largestExtent};
	nx = extentFactors[0]*(double(x)/double(voldimsGlobal[0]) - 0.5);
	ny = extentFactors[1]*(double(y)/double(voldimsGlobal[1]) - 0.5);
	nz = extentFactors[2]*(double(z)/double(voldimsGlobal[2]) - 0.5);
}

/* Converts `i` to a normalized coordinate (`nx`, `ny`, `nz`).
 * @i <BinaryVolume> index of interest
 * @nx the normalized *x* coordinate
 * @ny the normalized *y* coordinate
 * @nz the normalized *z* coordinate
 *
 * Converts `i` to a normalized coordinate (`nx`, `ny`, `nz`) according to `voldimsGlobal`.
 */
void normalizedCoord(voxelType i, double &nx, double &ny, double &nz){
	voxelType x(0), y(0), z(0);
	xyz(i, x, y, z, voldimsGlobal);
	normalizedCoord(x, y, z, nx, ny, nz);
}

string selectedString(){
	streamsize precDefault(4);
	double len(backbonesGlobal[ccSelected][bbSelected].length(voxdimsGlobal, voldimsGlobal)),
	vol(voxdimsGlobal[0]*voxdimsGlobal[1]*voxdimsGlobal[2]*volumesGlobal[ccSelected][bbSelected]);
	return "cc " + makeString(ccSelected) + ", bb " + makeString(bbSelected) + ", v " + makeString(vertSelected)
	+ ", vol=" + makeString(vol, precDefault) + "cu." + lengthUnitGlobal
	+ " (" + makeString(volumesGlobal[ccSelected][bbSelected]) + ")"
	+ " l=" + makeString(len, precDefault) + lengthUnitGlobal
	+ " <r_lv>=" + makeString(radFromVolLen(vol, len), precDefault) + lengthUnitGlobal
	+ " <r_so>=" + makeString(aveRadSolidGlobal[ccSelected][bbSelected], precDefault) + lengthUnitGlobal
	+ " <r_su>=" + makeString(aveRadSurfGlobal[ccSelected][bbSelected], precDefault) + lengthUnitGlobal;
}

/* Finds the closest <Backbone> to (`x`, `y`, `z`).
 * @x the normalized *x*-coordinate of the point of interest
 * @y the normalized *y*-coordinate of the point of interest
 * @z the normalized *z*-coordinate of the point of interest
 *
 * Finds the closest <Backbone> to (`x`, `y`, `z`), as long as the first global <Backbone> is not empty and the point is not too far away
 */
void findClosestSegment(double x, double y, double z){
	segIsSelected = false;
	if(backbonesGlobal.empty())
		return;
	if(backbonesGlobal[0].empty())
		return;
	if(backbonesGlobal[0][0].empty())
		return;
	if(2.0 < x || x < -2.0 || 2.0 < y || y < -2.0 || 2.0 < z || z < -2.0 || x != x || y != y || z != z)
		return;
	
	double nx(0.0), ny(0.0), nz(0.0);
	normalizedCoord(backbonesGlobal[0][0][0], nx, ny, nz);
	double sepSq(separationSq3D(x, y, z, nx, ny, nz));
	ccSelected = bbSelected = vertSelected = 0;
	for(unsigned int i(0); i < backbonesGlobal.size(); i++){
		for(unsigned int j(0); j < backbonesGlobal[i].size(); j++){
			for(unsigned int k(0); k < backbonesGlobal[i][j].size(); k++){
				normalizedCoord(backbonesGlobal[i][j][k], nx, ny, nz);
				double ss(separationSq3D(x, y, z, nx, ny, nz));
				if(sepSq > ss){
					sepSq = ss;
					ccSelected = i;
					bbSelected = j;
					vertSelected = k;
				}
			}
		}
	}
	segIsSelected = true;
}

/////////////// start wxMinSurfTests.h

#if defined(__WXMSW__) || defined(__WINDOWS__)
#include <windows.h>
#endif

#if defined(__WXMAC__)
//#	ifdef __DARWIN__
#		include <OpenGL/gl.h>
#		include <OpenGL/glu.h>
//#	else
//#		include <OpenGL/gl.h> // doh added directory
//#		include <OpenGL/glu.h> // doh added directory
//#	endif
#else
#	include <GL/gl.h>
#	include <GL/glu.h>
#endif

// the maximum number of vertices to accommodate
#define MAXVERTS	 1000000

// For compilers that support precompilation, includes "wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#if !wxUSE_GLCANVAS
#error "OpenGL required: set wxUSE_GLCANVAS to 1 and rebuild the library"
#endif

#include "wx/glcanvas.h"


future<void> fut, futSelected;

string selectedStatusString("");

GLboolean g_use_vertex_arrays = GL_FALSE; // Mac requires false
GLboolean g_doubleBuffer = GL_TRUE;

GLclampf clearColor[] = {1.0, 1.0, 1.0, 0.0};
GLdouble centerColor[] = {0.0, 1.0, 0.0, 1.0},
	selectColor[] = {0.0, 0.0, 1.0, 1.0},
	highlightColor1[4] = {1.0, 0.5, 0.5, 1.0},
	highlightColor2[4] = {0.5, 1.0, 0.5, 1.0},
	highlightColor3[4] = {0.5, 0.5, 1.0, 1.0},
	intrasphereColor[4] = {1.0, 0.0, 0.0, 0.2},
	intersphereColor[4] = {0.0, 1.0, 0.0, 0.2};

#include "FontBitmaps.h"

/* Writes the number `s` at (`x`, `y`, `z`) in color (`r`, `g`, `b`, `a`).
 * @s a string of a number that can include 0 through 9, and/or a negative sign, and/or a decimal point
 * @x the *x* coordinate to start drawing `s`
 * @y the *y* coordinate to start drawing `s`
 * @z the *z* coordinate to start drawing `s`
 * @r the red value of the color for drawing `s`
 * @g the green value of the color for drawing `s`
 * @b the blue value of the color for drawing `s`
 * @a the alpha value of the color for drawing `s`
 *
 * Writes the number converted to string `s` at position (`x`, `y`, `z`) with color variables (`r`, `g`, `b`, `a`).
 */
void writeNumber(string s, double x, double y, double z, double r = 0.0, double g = 0.0, double b = 0.0, double a = 1.0){
	glColor4d(r, g, b, a);
	glRasterPos3d(x, y, z);
	if(labelStatus%3 == 2){
		for(unsigned int i(0); i < s.length(); i++){
			switch(s[i]){
				case '-':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitDash16f); break;
				case '.':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitPeriod16f); break;
				case '0':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitZero16f); break;
				case '1':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitOne16f); break;
				case '2':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitTwo16f); break;
				case '3':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitThree16f); break;
				case '4':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitFour16f); break;
				case '5':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitFive16f); break;
				case '6':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitSix16f); break;
				case '7':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitSeven16f); break;
				case '8':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitEight16f); break;
				case '9':
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitNine16f); break;
				default:
					glBitmap(16, 16, 0.0f, 0.0f, 16.0f, 0.0f, bitSpace16f); break;
			}
		}
	}else if(labelStatus%3 == 0){
		for(unsigned int i(0); i < s.length(); i++){
			switch(s[i]){
				case '-':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitDash8f); break;
				case '.':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitPeriod8f); break;
				case '0':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitZero8f); break;
				case '1':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitOne8f); break;
				case '2':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitTwo8f); break;
				case '3':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitThree8f); break;
				case '4':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitFour8f); break;
				case '5':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitFive8f); break;
				case '6':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitSix8f); break;
				case '7':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitSeven8f); break;
				case '8':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitEight8f); break;
				case '9':
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitNine8f); break;
				default:
					glBitmap(8, 8, 0.0f, 0.0f, 8.0f, 0.0f, bitSpace8f); break;
			}
		}
	}
}

void writeNaturalNumber(unsigned int n, double x, double y, double z, double r = 0.0, double g = 0.0, double b = 0.0, double a = 1.0){
	writeNumber(makeString(n, 8), x, y, z, r, g, b, a);
}

void writeNumber(double w, double x, double y, double z, double r = 0.0, double g = 0.0, double b = 0.0, double a = 1.0){
	writeNumber(makeString(w, 8), x, y, z, r, g, b, a);
}

// normalize vector components in place
void normalize(GLfloat &x, GLfloat &y, GLfloat &z){
	GLfloat mag(sqrt(x*x + y*y + z*z));
	x /= mag;
	y /= mag;
	z /= mag;
}

class MyFrame;

// Define a new application type
class MyApp : public wxApp
{
public:
	virtual bool OnInit();// doh moved override: wxOVERRIDE;
	MyFrame *myFrame;
	
	//virtual void OnInitCmdLine(wxCmdLineParser& parser) wxOVERRIDE;
	//virtual bool OnCmdLineParsed(wxCmdLineParser& parser) wxOVERRIDE;
};

wxDECLARE_APP(MyApp);

// The OpenGL-enabled canvas
class TestGLCanvas : public wxGLCanvas{

private:
	wxGLContext* m_glRC;
	
	GLfloat m_verts[MAXVERTS][3],
		m_norms[MAXVERTS][3],
		m_cols[MAXVERTS][4];
	
	GLfloat m_numberverts[MAXVERTS][3];
	double m_numbers[MAXVERTS];
	GLfloat m_numbercols[MAXVERTS][4];
	
	GLint m_numverts,
		m_numlineverts,
		m_numnumbers;
	
	int winSize;
	
	GLfloat m_xrot;
	GLfloat m_yrot;
	
	GLfloat m_zoom,
		m_xtrans,
		m_ytrans,
		m_ztrans,
		m_eyex,
		m_eyey,
		m_eyez,
		m_centerx,
		m_centery,
		m_centerz,
		m_upx,
		m_upy,
		m_upz,
		m_rightx,
		m_righty,
		m_rightz,
		centerToEye[3], // for normals
		zNear,
		zFar;
	
	GLdouble ray[6];
	
	GLboolean modifierIsDown;
	
	void updateTranslation(double dx, double dy, double dz);
	void positionEye(double dx, double dy, double dz);
	void processClick(float x, float y);
	
	int getNewSize(int oldSize, int x, int y){
		if(x < y)
			return x;
		return y;
	}
	
	wxDECLARE_NO_COPY_CLASS(TestGLCanvas);
	wxDECLARE_EVENT_TABLE();
	
public:
	TestGLCanvas(wxWindow *parent,
				 wxWindowID id = wxID_ANY,
				 int *gl_attrib = NULL);
	
	virtual ~TestGLCanvas();
	
	void OnPaint(wxPaintEvent& WXUNUSED(event));
	
	void OnSize(wxSizeEvent& event){
		if (!IsShownOnScreen())
			return;
		
		SetCurrent(*m_glRC); // This is normally only necessary if there is more than one wxGLCanvas or more than one wxGLContext in the application.
	
		int newSize(getNewSize(winSize, event.GetSize().x, event.GetSize().y));
		//SetSize(newSize, newSize);
	
		// It's up to the application code to update the OpenGL viewport settings.
		// This is OK here only because there is only one canvas that uses the
		// context. See the cube sample for that case that multiple canvases are
		// made current with one context.
		// glViewport(0, 0, event.GetSize().x, event.GetSize().y);
		glViewport(0, 0, newSize, newSize);
		winSize = newSize;
	
		// doh wants a Refresh
		Refresh(false);
	}
	
	void OnChar(wxKeyEvent& event){
		switch(event.GetKeyCode()){
			//case WXK_ESCAPE:
			//	wxTheApp->ExitMainLoop();
			//	return;
			case WXK_LEFT:
				m_yrot -= 15.0;
				break;
			case WXK_RIGHT:
				m_yrot += 15.0;
				break;
			case WXK_UP:
				m_xrot += 15.0;
				break;
			case WXK_DOWN:
				m_xrot -= 15.0;
				break;
			case 'l': case 'L':
				if(modifierIsDown)
					labelStatus += 2;
				else
					labelStatus++;
				break;
			default:
				event.Skip();
				return;
		}
		Refresh(false);
	}
	
	void OnKeyDown(wxKeyEvent& event){
		switch(event.GetKeyCode()){
			case WXK_SHIFT:
				modifierIsDown = true;
				Refresh(false);
				return;
			default:
				event.Skip();
				return;
		}
	}
	
	void OnKeyUp(wxKeyEvent& event){
		switch(event.GetKeyCode()){
			case WXK_SHIFT:
				modifierIsDown = false;
				Refresh(false);
				return;
			default:
				event.Skip();
				return;
		}
	}
	
	void OnMouseEvent(wxMouseEvent& event){
		static int dragging = 0;
		static float last_x, last_y, first_x, first_y;
	
		// Allow default processing to happen, or else the canvas cannot gain focus
		// (for key events).
		event.Skip();
	
		if(event.LeftIsDown()){
			if (!dragging){
				first_x = event.GetX();
				first_y = event.GetY();
				dragging = 1;
			}
			else{
				//updateTranslation(event.GetX() - last_x, event.GetY() - last_y, 0.0);
				positionEye(event.GetX() - last_x, event.GetY() - last_y, 0.0);
			
			}
			last_x = event.GetX();
			last_y = event.GetY();
		}else{
			if(wxGetApp().myFrame != NULL && first_x == event.GetX() && first_y == event.GetY())
				processClick(first_x, first_y);
			dragging = 0;
		}
	
		if(event.GetWheelRotation() != 0)
			positionEye(0.0, 0.0, 1.0*event.GetWheelRotation()/event.GetWheelDelta());
	
		Refresh(false);
	}
	
	// surface is ordered (x0, y0, z0, nx0, ny0, nz0, x1,..., nzi) for i triangles
	void newSurface(const std::vector<double> &surface);
	
	// scaffold is ordered (x0, y0, z0, r0, g0, b0, a0, x1,..., ai) for i vertices for GL_LINE_STRIP
	void newScaffold(const std::vector<double> &surface);
	
	// points is ordered (x0, y0, z0, r0, g0, b0, a0, x1,..., ai) for i vertices for GL_POINTS
	void newPoints(const std::vector<double> &points);
	
	// lines is ordered (x0, y0, z0, r0, g0, b0, a0, x1,..., ai) for i vertices for GL_POINTS
	void newLines(const std::vector<double> &lines){
		m_numlineverts = 0;
		for(unsigned int i(0); i < lines.size()/7 && m_numlineverts <= MAXVERTS; i++, m_numlineverts++){
			for(unsigned int j(0); j < 3; j++){
				m_verts[m_numverts + m_numlineverts][j] = lines[7*i + j];
				m_cols[m_numverts + m_numlineverts][j] = lines[7*i + j + 3];
			}
			m_cols[m_numverts + m_numlineverts][3] = lines[7*i + 6];
		}
		Refresh(false);
	}
	
	// adds line with color; expects (x, y, z, r, g, b, a)
	void addLine(const std::vector<double> &line){
		unsigned int lineStart(m_numverts + m_numlineverts);
		if(line.empty() || m_numlineverts + line.size()/7 + 2 > MAXVERTS)
			return;
		
		m_numlineverts += line.size()/7 + 2;
		for(unsigned int j(0); j < 3; j++){ // blank transparent tip
			m_verts[lineStart][j] = line[j];
			m_cols[lineStart][j] = 0;
		}
		m_cols[lineStart][3] = 0;
		for(unsigned int i(0); i < line.size()/7 && m_numlineverts <= MAXVERTS; i++){
			for(unsigned int j(0); j < 3; j++){
				m_verts[lineStart + 1 + i][j] = line[7*i + j];
				m_cols[lineStart + 1 + i][j] = line[7*i + j + 3];
			}
			m_cols[lineStart + 1 + i][3] = line[7*i + 6];
		}
		for(unsigned int j(0); j < 3; j++){ // blank transparent tip
			m_verts[lineStart + line.size()/7 + 1][j] = line[line.size() - 7 + j];
			m_cols[lineStart + line.size()/7 + 1][j] = 0;
		}
		m_cols[lineStart + line.size()/7 + 1][3] = 0;
	}
	
	// update colors according to uc
	void updateCols(const map<unsigned int, std::vector<double> > &uc){
		for(map<unsigned int, vector<double> >::const_iterator it(uc.begin()); it != uc.end(); it++){
			unsigned int i(it->first);
			for(short j(0); j < 4; j++)
				m_cols[i][j] = it->second[j];
		}
		Refresh(false);
	}
	
	// update all colors to (r, g, b, a)
	void updateCols(double r, double g, double b, double a){
		for(unsigned int i(0); i < m_numverts; i++){
			m_cols[i][0] = r;
			m_cols[i][1] = g;
			m_cols[i][2] = b;
			m_cols[i][3] = a;
		}
		Refresh(false);
	}
	
	// vector `numbers` in groups of 8: (x, y, z) (r, g, b, a) <number>
	void newNumbers(const std::vector<double> &numbers){
		m_numnumbers = 0;
		for(unsigned int i(0); i < numbers.size()/8 && m_numnumbers <= MAXVERTS; i++, m_numnumbers++){
			for(unsigned int j(0); j < 3; j++){
				m_numberverts[m_numnumbers][j] = numbers[8*i + j];
				m_numbercols[m_numnumbers][j] = numbers[8*i + j + 3];
			}
			normalize(m_norms[m_numnumbers][0], m_norms[m_numnumbers][1], m_norms[m_numnumbers][2]);
			m_numbercols[m_numnumbers][3] = numbers[8*i + 6];
			m_numbers[m_numnumbers] = numbers[8*i + 7];
		}
		Refresh(false);
	}
	
	// add the number w and (x, y, z) with color (r, g, b, a) [default is black]
	void addNumber(double w, double x, double y, double z, double r = 0.0, double g = 0.0, double b = 0.0, double a = 1.0){
		GLdouble numbers[] = {x, y, z, r, g, b, a, w};
		for(unsigned int j(0); j < 3; j++){
			m_numberverts[m_numnumbers][j] = numbers[j];
			m_numbercols[m_numnumbers][j] = numbers[j + 3];
		}
		normalize(m_norms[m_numnumbers][0], m_norms[m_numnumbers][1], m_norms[m_numnumbers][2]);
		m_numbercols[m_numnumbers][3] = numbers[6];
		m_numbers[m_numnumbers] = numbers[7];
		m_numnumbers++;
		Refresh(false);
	}
	
	//void LoadSurface(const wxString& filename);
	
	// initialize material properties
	void InitMaterials(){
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(5);
		glPointSize(4);
	}
	
	// initialize GL properties
	void InitGL(){
		// Make the new context current (activate it for use) with this canvas.
		SetCurrent(*m_glRC);
	
		glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);
		glShadeModel(GL_SMOOTH);
		glEnable(GL_DEPTH_TEST);
		InitMaterials();
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-1.0, 1.0, -1.0, 1.0, zNear, zFar);
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef( 0.0, 0.0, -zNear - 1.0 );
		
		if (g_use_vertex_arrays){
			glVertexPointer(3, GL_FLOAT, 0, m_verts);
			glNormalPointer(GL_FLOAT, 0, m_norms);
			glColorPointer(4, GL_FLOAT, 0, m_cols);
		
			glVertexPointer(3, GL_FLOAT, m_numverts, m_verts);
			glColorPointer(4, GL_FLOAT, m_numverts, m_cols);
		
			glEnable(GL_VERTEX_ARRAY);
			glEnable(GL_NORMAL_ARRAY);
			glEnable(GL_COLOR_ARRAY);
		}
		
		glEnable(GL_FOG);
		//glFogf(GL_FOG_DENSITY, 0.5f);
		glFogfv(GL_FOG_COLOR, clearColor);
		glFogi(GL_FOG_MODE, GL_LINEAR);
		positionEye(0.0, 0.0, 0.0);
		InitMaterials();
	}
};

// The frame containing the GL canvas
class MyFrame : public wxFrame
{
public:
	MyFrame(wxFrame *frame,
			const wxString& title,
			const wxPoint& pos = wxDefaultPosition,
			const wxSize& size = wxDefaultSize,
			long style = wxDEFAULT_FRAME_STYLE);
	
	virtual ~MyFrame();
	
	TestGLCanvas *m_canvas;
	
	// updates status message with `s`, keeping `selectedStatusString`
	void updateStatus(std::string s){
		wxGetApp().myFrame->SetStatusText(wxString(s) + " " + selectedStatusString);
	}
	
	private :
	// Intercept menu commands
	void OnExit( wxCommandEvent& WXUNUSED(event) ){
		Close(true);// true is to force the frame to close
	}
	
	wxDECLARE_EVENT_TABLE();
};

///////////////////////// end wxMinSurfTests.h

///////////////////////// start wxMinSurfTests.cpp

// adapted from isosurf.cpp by Brian Paul and Wolfram Gloger

// For compilers that support precompilation, includes "wx.h".

#include "wx/timer.h"
#include "wx/glcanvas.h"
#include "wx/math.h"
#include "wx/log.h"
#include "wx/cmdline.h"
#include "wx/wfstream.h"
#include "wx/zstream.h"
#include "wx/txtstrm.h"

//#include "wxMinSurfTests.h"

//---------------------------------------------------------------------------
// MyApp
//---------------------------------------------------------------------------

wxIMPLEMENT_APP(MyApp);

void sphereCoarsenTest();

bool MyApp::OnInit(){
	if ( !wxApp::OnInit() )
		return false;
	
	initBitmapNumbers();
	
#if TRY_WX==1
	
	
	// Create the main frame window // Adapted from wxWidgets OpenGL Isosurf Sample...
	myFrame = new MyFrame(NULL, wxT("Angicart++"), wxPoint(50, 50), wxSize(640, 710)); // taller for menu and status
	
	
	fut = std::async(sphereCoarsenTest);// wxTest sphereCoarsenTest
	
	return true;
#endif
	return false;
}

wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
EVT_MENU(wxID_EXIT, MyFrame::OnExit)
wxEND_EVENT_TABLE()

MyFrame::MyFrame(wxFrame *frame, const wxString& title, const wxPoint& pos, const wxSize& size, long style)
: wxFrame(frame, wxID_ANY, title, pos, size, style), m_canvas(NULL){
	
	wxMenu *fileMenu = new wxMenu; //menubar
	
	fileMenu->Append(wxID_EXIT, wxT("E&xit"));
	wxMenuBar *menuBar = new wxMenuBar;
	menuBar->Append(fileMenu, wxT("&File"));
	SetMenuBar(menuBar);
	
#ifdef __WXMSW__
	int *gl_attrib = NULL;
#else
	int gl_attrib[20] = { WX_GL_RGBA, WX_GL_MIN_RED, 1, WX_GL_MIN_GREEN, 1,
		WX_GL_MIN_BLUE, 1, WX_GL_DEPTH_SIZE, 1,
		WX_GL_DOUBLEBUFFER,
#	if defined(__WXMAC__) || defined(__WXQT__)
		GL_NONE };
#	else
	None };
#	endif
#endif

	if (!g_doubleBuffer){
		wxLogWarning("Disabling double buffering");
#ifdef __WXGTK__
		gl_attrib[9] = None;
#endif
		g_doubleBuffer = GL_FALSE;
	}

	m_canvas = new TestGLCanvas(this, wxID_ANY, gl_attrib);

	Show(true); // show frame
	Raise();

	m_canvas->InitGL();

	CreateStatusBar();
}

MyFrame::~MyFrame(){ delete m_canvas; }

//---------------------------------------------------------------------------
// TestGLCanvas
//---------------------------------------------------------------------------

wxBEGIN_EVENT_TABLE(TestGLCanvas, wxGLCanvas)
EVT_SIZE(TestGLCanvas::OnSize)
EVT_PAINT(TestGLCanvas::OnPaint)
EVT_CHAR(TestGLCanvas::OnChar)
EVT_KEY_DOWN(TestGLCanvas::OnKeyDown)
EVT_KEY_UP(TestGLCanvas::OnKeyUp)
EVT_MOUSE_EVENTS(TestGLCanvas::OnMouseEvent)
wxEND_EVENT_TABLE()

TestGLCanvas::TestGLCanvas(wxWindow *parent, wxWindowID id, int* gl_attrib)
: wxGLCanvas(parent, id, gl_attrib){
	m_xrot = 0;
	m_yrot = 0;
	m_xtrans = 0;
	m_ytrans = 0;
	m_ztrans = 0.0f;
	m_zoom = -0.3f;
	
	m_eyex = 0.0f;
	m_eyey = 0.0f;
	m_eyez = 5.0f;
	m_centerx = m_centery = m_centerz = 0.0f;
	m_upx = 0.0;
	m_upy = 1.0;
	m_upz = 0.0;
	m_rightx = 1.0;
	m_righty = 0.0;
	m_rightz = 0.0;
	
	zNear = 5.0f;
	zFar = 6.0f;
	
	GetSize(&winSize, &winSize);
	SetSize(winSize, winSize);
	
	for(unsigned int i(0); i < 6; i++)
		ray[i] = 0.0f;
	
	m_numverts = 0;
	
	// Explicitly create a new rendering context instance for this canvas.
	m_glRC = new wxGLContext(this);
}

TestGLCanvas::~TestGLCanvas(){ delete m_glRC; }

// re-initializes surface
void TestGLCanvas::newSurface(const std::vector<double> &surface){
	m_numverts = 0;
	for(unsigned int i(0); i < surface.size()/6 && m_numverts <= MAXVERTS; i++){
		unsigned int j(6*i);
		for(unsigned int k(0); k < 3; k++){
			m_verts[m_numverts][k] = surface[j + k];
			m_norms[m_numverts][k] = surface[j + 3 + k];
		}
		m_cols[m_numverts][0] = GLfloat(i)/GLfloat(surface.size()/6);
		m_cols[m_numverts][1] = GLfloat(surface.size()/6 - i)/GLfloat(surface.size()/6);
		m_cols[m_numverts][2] = GLfloat(i)/GLfloat(surface.size()/6);
		m_cols[m_numverts][3] = 0.5;
		m_numverts++;
	}
	Refresh(false);
}

// re-initializes scaffold
void TestGLCanvas::newScaffold(const std::vector<double> &scaffold){
	m_numverts = 0;
	for(unsigned int i(0); i < scaffold.size()/7 && m_numverts <= MAXVERTS; i++, m_numverts++){
		for(unsigned int j(0); j < 3; j++){
			m_verts[m_numverts][j] = scaffold[7*i + j];
			m_cols[m_numverts][j] = scaffold[7*i + j + 3];
		}
		m_cols[m_numverts][3] = scaffold[7*i + 6];
	}
	Refresh(false);
}

// re-initializes points
void TestGLCanvas::newPoints(const std::vector<double> &points){
	unsigned int newNumverts((unsigned int)points.size()/7);
	if(newNumverts > m_numverts){ // move line verts to lower indices (start at back)
		for(unsigned int k(0); k < m_numlineverts; k++){
			unsigned int i(m_numlineverts - k - 1);
			for(short j(0); j < 3; j++){
				m_verts[newNumverts + i][j] = m_verts[m_numverts + i][j];
				m_norms[newNumverts + i][j] = m_norms[m_numverts + i][j];
				m_cols[newNumverts + i][j] = m_cols[m_numverts + i][j];
			}
			m_cols[newNumverts + i][3] = m_cols[m_numverts + i][3];
		}
	}else if(newNumverts < m_numverts){ // move line verts to higher indices (start at front)
		for(unsigned int i(0); i < m_numlineverts; i++){
			for(short j(0); j < 3; j++){
				m_verts[newNumverts + i][j] = m_verts[m_numverts + i][j];
				m_norms[newNumverts + i][j] = m_norms[m_numverts + i][j];
				m_cols[newNumverts + i][j] = m_cols[m_numverts + i][j];
			}
			m_cols[newNumverts + i][3] = m_cols[m_numverts + i][3];
		}
	}
	
	// update points
	m_numverts = 0;
	for(unsigned int i(0); i < points.size()/7 && m_numverts <= MAXVERTS; i++, m_numverts++){
		for(unsigned int j(0); j < 3; j++){
			m_verts[m_numverts][j] = points[7*i + j];
			m_norms[m_numverts][j] = uniformRN() - 0.5;
			m_cols[m_numverts][j] = points[7*i + j + 3];
		}
		normalize(m_norms[m_numverts][0], m_norms[m_numverts][1], m_norms[m_numverts][2]);
		m_cols[m_numverts][3] = points[7*i + 6];
	}
	Refresh(false);
}


void TestGLCanvas::OnPaint(wxPaintEvent& WXUNUSED(event)){
	wxPaintDC dc(this); // This is a dummy, to avoid an endless succession of paint messages.  OnPaint handlers must always create a wxPaintDC.

	SetCurrent(*m_glRC); // This is normally only necessary if there is more than one wxGLCanvas or more than one wxGLContext in the application.

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	
	// draw
	if (g_use_vertex_arrays){
		glDrawArrays(GL_POINTS, 0, m_numverts);
		glDrawArrays(GL_LINE_STRIP, m_numverts, m_numlineverts);
	}else{
		glBegin(GL_POINTS);
		for(int i(0); i < m_numverts; i++){
			glColor4fv(m_cols[i]);
			glNormal3fv(m_norms[i]);
			glVertex3fv(m_verts[i]);
		}
		glEnd();
		glBegin(GL_LINE_STRIP);
		for(int i(m_numverts); i < m_numverts + m_numlineverts; i++){
			glColor4fv(m_cols[i]);
			glVertex3fv(m_verts[i]);
		}
		glEnd();
	}

	// point at focus
	glBegin(GL_POINTS);
	glColor4dv(centerColor);
	glVertex3d(m_centerx, m_centery, m_centerz);
	glEnd();

	// trace ray
	if(ray[0] == ray[3] || ray[1] == ray[4] || ray[2] == ray[5]){
		glBegin(GL_POINTS);
		glColor4dv(selectColor);
		glVertex3dv(&ray[0]);
		glEnd();
	}else{
		glBegin(GL_LINES);
		glColor4dv(selectColor);
		glVertex3dv(&ray[0]);
		glColor4dv(selectColor);
		glVertex3dv(&ray[3]);
		glEnd();
		glBegin(GL_POINTS);
		glColor4dv(selectColor);
		glVertex3d((ray[3] + ray[0])/2.0, (ray[4] + ray[1])/2.0, (ray[5] + ray[2])/2.0);
		glEnd();
	}

	// selected segment
	if(segIsSelected){
		glBegin(GL_POINTS);
		for(unsigned int i(0); i < backbonesGlobal[ccSelected][bbSelected].size(); i++){
			double nx(0.0), ny(0.0), nz(0.0);
			normalizedCoord(backbonesGlobal[ccSelected][bbSelected][i], nx, ny, nz);
			double offsetInc(0.01);
			if(i%3 == 0){
				glColor4dv(highlightColor1);
				glVertex3d(nx + offsetInc, ny, nz);
				glColor4dv(highlightColor1);
				glVertex3d(nx - offsetInc, ny, nz);
			}
			else if(i%3 == 1){
				glColor4dv(highlightColor2);
				glVertex3d(nx, ny + offsetInc, nz);
				glColor4dv(highlightColor2);
				glVertex3d(nx, ny - offsetInc, nz);
			}else{
				glColor4dv(highlightColor3);
				glVertex3d(nx, ny, nz + offsetInc);
				glColor4dv(highlightColor3);
				glVertex3d(nx, ny, nz - offsetInc);
			}
			glVertex3d(nx, ny, nz);
		}
		glEnd();
	}
	if(updateSelectedStatus){
		updateSelectedStatus = false;
		if(segIsSelected){
			selectedStatusString = selectedString();
			wxGetApp().myFrame->updateStatus("Selected");
		}else
			wxGetApp().myFrame->updateStatus("");
	}

	GLfloat position[] = {m_eyex, m_eyey, m_eyez, 1.0f};
	glLightfv(GL_LIGHT0, GL_POSITION, position);

	// numbers
	for(int i(0); i < m_numnumbers; i++){
		writeNumber(m_numbers[i], m_numberverts[i][0], m_numberverts[i][1], m_numberverts[i][2],
					m_numbercols[i][0], m_numbercols[i][1], m_numbercols[i][2], m_numbercols[i][3]);
	}

	glPopMatrix();
	SwapBuffers();
}





// cross product: a x b = c
inline void crossProduct(const GLfloat &ax, const GLfloat &ay, const GLfloat &az, const GLfloat &bx, const GLfloat &by, const GLfloat &bz, GLfloat &cx, GLfloat &cy, GLfloat &cz){
	cx = ay*bz - az*by;
	cy = az*bx - ax*bz;
	cz = ax*by - ay*bx;
}

// a is eye, b is center, c is vector from eye to center
inline void ahead(const GLfloat &ax, const GLfloat &ay, const GLfloat &az, const GLfloat &bx, const GLfloat &by, const GLfloat &bz, GLfloat &cx, GLfloat &cy, GLfloat &cz){
	cx = bx - ax;
	cy = by - ay;
	cz = bz - az;
}

// magnitude of triplet
inline GLfloat magnitude(GLfloat x, GLfloat y, GLfloat z){
	return sqrt(x*x + y*y + z*z);
};

// moves the camera or center
void TestGLCanvas::positionEye(double dx, double dy, double dz){
	
	GLfloat aheadx(0), aheady(0), aheadz(0);
	ahead(m_eyex, m_eyey, m_eyez, m_centerx, m_centery, m_centerz, aheadx, aheady, aheadz);
	GLfloat r(magnitude(aheadx, aheady, aheadz));
	
	if(modifierIsDown){ // move center
		// approximate movement in tangent plane
		GLfloat scaleCoefx(-0.003f), scaleCoefy(-0.003f), scaleCoefz(-0.003f);
		m_centerx += scaleCoefx*m_rightx*dx - scaleCoefy*m_upx*dy;
		m_centery += scaleCoefx*m_righty*dx - scaleCoefy*m_upy*dy;
		m_centerz += scaleCoefx*m_rightz*dx - scaleCoefy*m_upz*dy;
		
		// adjust separation of eye and center to make a rotation
		ahead(m_eyex, m_eyey, m_eyez, m_centerx, m_centery, m_centerz, aheadx, aheady, aheadz);
		GLfloat rNew(magnitude(aheadx, aheady, aheadz));
		m_centerx = m_eyex + r*aheadx/rNew;
		m_centery = m_eyey + r*aheady/rNew;
		m_centerz = m_eyez + r*aheadz/rNew;
		
		// move in direction of ahead (could be moved before updating dx and dy movement)
		normalize(aheadx, aheady, aheadz);
		m_centerx += scaleCoefz*aheadx*dz;
		m_centery += scaleCoefz*aheady*dz;
		m_centerz += scaleCoefz*aheadz*dz;
		m_eyex += scaleCoefz*aheadx*dz;
		m_eyey += scaleCoefz*aheady*dz;
		m_eyez += scaleCoefz*aheadz*dz;
	}else{ // move eye
		// approximate movement in tangent plane
		GLfloat scaleCoefx(0.1f), scaleCoefy(0.1f), scaleCoefz(0.1f);
		m_eyex -= scaleCoefx*m_rightx*dx - scaleCoefy*m_upx*dy;
		m_eyey -= scaleCoefx*m_righty*dx - scaleCoefy*m_upy*dy;
		m_eyez -= scaleCoefx*m_rightz*dx - scaleCoefy*m_upz*dy;
		ahead(m_eyex, m_eyey, m_eyez, m_centerx, m_centery, m_centerz, aheadx, aheady, aheadz);
		
		// adjust separation of eye and center to make a rotation
		GLfloat rNew(magnitude(m_centerx - m_eyex, m_centery - m_eyey, m_centerz - m_eyez));
		m_eyex = m_centerx - r*aheadx/rNew;
		m_eyey = m_centery - r*aheady/rNew;
		m_eyez = m_centerz - r*aheadz/rNew;
		
		// zoom eye
		m_zoom -= scaleCoefz*dz;
		double zoomFact(pow(10.0, m_zoom));
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		zNear = r - sqrt(3.0);
		zFar = r + sqrt(3.0);
		glFrustum(-zoomFact, zoomFact, -zoomFact, zoomFact, zNear, zFar);
		glFogf(GL_FOG_START, r - sqrt(3.0));
		glFogf(GL_FOG_END, r + sqrt(3.0));
	}
	
	// adjust up to be perpendicular to ahead
	ahead(m_eyex, m_eyey, m_eyez, m_centerx, m_centery, m_centerz, aheadx, aheady, aheadz);
	normalize(aheadx, aheady, aheadz);
	m_upx -= m_upx*aheadx;
	m_upy -= m_upy*aheady;
	m_upz -= m_upz*aheadz;
	normalize(m_upx, m_upy, m_upz);
	
	// update MODELVIEW
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(m_eyex, m_eyey, m_eyez, m_centerx, m_centery, m_centerz, m_upx, m_upy, m_upz);
	
	// adjust right (rightward direction) to be perpendicular to ahead and up
	crossProduct(aheadx, aheady, aheadz, m_upx, m_upy, m_upz, m_rightx, m_righty, m_rightz);
	normalize(m_rightx, m_righty, m_rightz);
}














void TestGLCanvas::updateTranslation(double dx, double dy, double dz){
	double radPerDeg(acos(-1.0)/180.0);
	double DX(dx*cos(m_yrot*radPerDeg) + dy*sin(m_xrot*radPerDeg)),
	DY(dy*cos(m_xrot*radPerDeg) + dx*sin(m_yrot*radPerDeg)),
	DZ(dz);
	m_xtrans += 0.01*DX;
	m_ytrans -= 0.01*DY;
	m_ztrans += 0.01*DZ;
}

void TestGLCanvas::processClick(float x, float y){
	int w, h;
	GetSize(&w, &h);
	y = h - y;
	GLdouble model[16], proj[16];	
	GLint viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, model);
	glGetDoublev(GL_PROJECTION_MATRIX, proj);
	glGetIntegerv(GL_VIEWPORT, viewport);
	if(ray[0] == ray[3] && ray[1] == ray[4] && ray[2] == ray[5]){
		segIsSelected = false;
		selectedStatusString = "";
		GLint suc = gluUnProject(x, y, 0.0f, model, proj, viewport, &ray[0], &ray[1], &ray[2]);
		suc = suc && gluUnProject(x, y, 1.0f, model, proj, viewport, &ray[3], &ray[4], &ray[5]);
	}else{
		GLdouble x0(0.0), y0(0.0), z0(0.0), x1(0.0), y1(0.0), z1(0.0);
		GLint suc = gluUnProject(x, y, -1.0f, model, proj, viewport, &x0, &y0, &z0);
		suc = suc && gluUnProject(x, y, 1.0f, model, proj, viewport, &x1, &y1, &z1);
		GLfloat m[3], n[3], b[3], c[3];
		for(unsigned int i(0); i < 3; i++){
			m[i] = ray[i + 3] - ray[i];
			b[i] = ray[i];
		}
		n[0] = x1 - x0;
		n[1] = y1 - y0;
		n[2] = z1 - z0;
		c[0] = x0;
		c[1] = y0;
		c[2] = z0;
		normalize(m[0], m[1], m[2]);
		normalize(n[0], n[1], n[2]);
		float mn(0.0), dm(0.0), ys(0.0);
		for(unsigned int i(0); i < 3; i++)
			mn += m[i]*n[i];
		if(abs(mn) >= 1.0) // parallel
			return;
		for(unsigned int i(0); i < 3; i++){
			dm += (b[i] - c[i])*m[i];
			ys += (b[i] - c[i])*(mn*m[i] - n[i]);
		}
		ys /= mn*mn - 1.0;
		for(unsigned int i(0); i < 3; i++)
			ray[i] = ray[i + 3] = (m[i]*(mn*ys - dm) + n[i]*ys + b[i] + c[i])/2.0;
		futSelected = std::async(findClosestSegment, ray[0], ray[1], ray[2]);
	}
	updateSelectedStatus = true;
	Refresh(false);
}



//////////////////////// end wxMinSurfTests.cpp


/* Finds the lowest unsigned integer that is not in `v`
 * @v the list of unsigned integers that should not be returned
 *
 * Starts at zero and increments by one until an unsigned integer is found that is not in `v`.
 * @return the lowest unsigned integer that is not in `v`
 */
unsigned int lowestNotIn(const vector<unsigned int> &v){
	unsigned int i(0);
	while(isIn(i, v))
		i++;
	return i;
}

/* Assigns a color to each <Backbone> in `backbones` that tries to be different from its neighbors.
 * @backbones the list of <Backbone> objects for color assignment
 * @branchpoints sets of adjacent backbones after segmentation (i.e., the  connectivity information)
 * @numCol number of different colors to use
 *
 * Assigns a color to each <Backbone> in `backbones` that tries to be different from its neighbors.  Iterates through `numCol` colors and attempts to resolve conflicts by looking for a color that is not yet in the neighborhood.
 * @return list of colors corresponding to each <Backbone> ins `backbones`
 */
vector<vector<vector<unsigned char> > > segColors(const vector<vector<Backbone<> > > &backbones,
		const vector<vector<vector<branchType> > > &branchpoints, unsigned int numCol = 9){
	vector<vector<vector<unsigned char> > > sc(branchpoints.size(), vector<vector<unsigned char> >());
	if(branchpoints.size() == 0){
		sc.push_back(vector<vector<unsigned char> >(1, vector<unsigned char>(4, 255)));
		rainbowColor(0, numCol, sc[0][0][0], sc[0][0][1], sc[0][0][2]);
		return sc;
	}
	
	if(numCol%2 == 0)
		numCol++;
	vector<vector<unsigned char> > cols(numCol, vector<unsigned char>(4, 255));
	for(unsigned int i(0); i < numCol; i++)
		rainbowColor(i, numCol, cols[i][0], cols[i][1], cols[i][2]);

	vector<map<voxelType, vector<branchType> > > adj(branchpoints.size(), map<voxelType, vector<branchType> >());
	vector<vector<unsigned int> > label(branchpoints.size(), vector<unsigned int>());
	for(unsigned int i(0); i < branchpoints.size(); i++){ // each cc
		// construct adjacency
		for(unsigned int j(0); j < branchpoints[i].size(); j++){ // each branchpoint in cc
			for(unsigned int k(0); k < branchpoints[i][j].size(); k++){ // each backbone in branchpoint in cc
				for(unsigned int m(0); m < k; m++)
					pushUnique(branchpoints[i][j][m], adj[i][branchpoints[i][j][k]]);
				for(unsigned int m(k + 1); m < branchpoints[i][j].size(); m++)
					pushUnique(branchpoints[i][j][m], adj[i][branchpoints[i][j][k]]);
			}
		}

		// assign labels
		label[i] = vector<unsigned int>(backbones[i].size(), 0);
		for(unsigned int j(0); j < backbones[i].size(); j++)
			label[i][j] = j%numCol;
		for(unsigned int j(0); j < backbones[i].size(); j++){
			vector<unsigned int> nhCols;
			for(unsigned int k(0); k < adj[i][j].size(); k++)
				pushUnique(label[i][adj[i][j][k]], nhCols);
			if(isIn(label[i][j], nhCols))
				label[i][j] = lowestNotIn(nhCols);
			if(label[i][j] >= numCol)
				label[i][j] = j%numCol; // revert to some color, even if it is a duplicate
		}
	}

	for(unsigned int i(0); i < branchpoints.size(); i++){
		sc[i] = vector<vector<unsigned char> >(backbones[i].size(), vector<unsigned char>(4, 255));
		for(unsigned int j(0); j < backbones[i].size(); j++){
			for(unsigned int k(0); k < 3; k++)
				sc[i][j][k] = cols[label[i][j]][k];
		}
	}

	return sc;
}

/* Returns a vector for the <Backbone>'s color, formatted for adding to wx
 * @backbone <Backbone> with vertebrae to define line
 * @col 4-vector to define line color
 *
 * Returns a vector for the <Backbone>'s color, formatted for adding to wx
 * @ return a vector for the <Backbone>'s color, formatted for adding to wx
 */
vector<double> lineCols(const Backbone<> &backbone, const vector<unsigned char> &col){
	vector<double> bc, fourTrans(4, 0);
	double nx, ny, nz;
	if(backbone.size() > 0){ // blank and transparent buffer point at backbone start
		normalizedCoord(backbone[0], nx, ny, nz);
		bc.push_back(nx);
		bc.push_back(ny);
		bc.push_back(nz);
		bc.insert(bc.end(), fourTrans.begin(), fourTrans.end());
	}
	for(unsigned int k(0); k < backbone.size(); k++){ // each vertebra
		normalizedCoord(backbone[k], nx, ny, nz);
		bc.push_back(nx);
		bc.push_back(ny);
		bc.push_back(nz);
		bc.push_back(col[0]/255.0);
		bc.push_back(col[1]/255.0);
		bc.push_back(col[2]/255.0);
		bc.push_back(col[3]/255.0);
	}
	if(backbone.size() > 0){ // blank and transparent buffer point at backbone end
		normalizedCoord(backbone.back(), nx, ny, nz);
		bc.push_back(nx);
		bc.push_back(ny);
		bc.push_back(nz);
		bc.insert(bc.end(), fourTrans.begin(), fourTrans.end());
	}
	return bc;
}

/* Constructs a vector to visualize all points in <BinaryVolume> `B`.
 * @B <BinaryVolume> that indicates points to visualize
 * @pm the map from <BinaryVolume> indices to index of visualized point
 * @sets sets of points that will be assigned particular colors
 * @setCols list of color properties for the sets: a single vector (r0, g0, b0, a0, r1, ..., aN)
 *
 * Constructs a vector to visualize all points in <BinaryVolume> `B`, with the specific points in `sets` assigned custom colors.
 * @return a vector to visualize all points in <BinaryVolume> `B`
 */
vector<double> pointsWithColors(const BinaryVolume &B, map<voxelType, unsigned int> &pm, const vector<vector<voxelType> > &sets,
		const vector<double> &setCols){
	map<voxelType, unsigned int> pointSets;
	for(voxelType i(B.findFirstAtOrAfter(0)); i < B.totalSize(); i = B.findFirstAtOrAfter(i + 1))
		pointSets[i] = (unsigned int)sets.size(); // default is not assigned to any set
	for(unsigned int i(0); i < sets.size(); i++){
		for(unsigned int j(0); j < sets[i].size(); j++)
			pointSets[sets[i][j]] = i; // assign label for those points that are in a set
	}
	vector<double> v(7*pointSets.size(), 1.0);
	unsigned int i(0);
	for(map<voxelType, unsigned int>::iterator it(pointSets.begin()); it != pointSets.end(); it++){
		normalizedCoord(it->first, v[7*i], v[7*i + 1], v[7*i + 2]);
		for(unsigned int j(3); j < 7; j++)
			v[7*i + j] = defaultPointColor[j - 3]; // set default color
		pm[it->first] = i;
		if(it->second == sets.size()) // has not been assigned to any set
			v[7*i + 6] = 0.2;
		else if(4*sets.size() == setCols.size()){ // set to the specially prescribed color
			v[7*i + 3] = setCols[4*it->second];
			v[7*i + 4] = setCols[4*it->second + 1];
			v[7*i + 5] = setCols[4*it->second + 2];
			v[7*i + 6] = setCols[4*it->second + 3];
		}else{ // default coloring
			v[7*i + 3] = 0.2;
			v[7*i + 4] = 0.6;
			v[7*i + 5] = 0.4;
			v[7*i + 6] = 0.4;
		}
		i++;
	}
	
	return v;
}

/* Constructs a vector to visualize all points in <BinaryVolume> `B`.
 * @B <BinaryVolume> that indicates points to visualize
 * @pm the map from <BinaryVolume> indices to index of visualized point
 *
 * Constructs a vector to visualize all points in <BinaryVolume> `B`, with the specific points in `sets` assigned custom colors.
 * @return a vector to visualize all points in <BinaryVolume> `B`
 */
vector<double> pointsWithColors(const BinaryVolume &B, map<voxelType, unsigned int> &pm){
	return pointsWithColors(B, pm, vector<vector<voxelType> >(), vector<double>());
}

/* Sets the points (with <BinaryVolume> indices in `v`) to the default color.
 * @v the <BinaryVolume> indices of the voxels of interest
 * @pm the map from <BinaryVolume> indices to index of visualized point
 *
 * Updates the colors of all the points specified in `v` by their <BinaryVolume> index to the default color.
 */
template <class T>
void setDefaultColor(const vector<T> &v, const map<voxelType, unsigned int> &pm){
	map<unsigned int, vector<double> > uc;
	for(unsigned i(0); i < v.size(); i++){
		map<voxelType, unsigned int>::const_iterator pmit(pm.find(v[i]));
		if(pmit == pm.end())
			continue;
		
		unsigned int pmIndex(pmit->second);
		for(unsigned int j(0); j < 4; j++)
			uc[pmIndex].push_back(defaultPointColor[j]);
	}
	wxGetApp().myFrame->m_canvas->updateCols(uc);
}

/* Sets the points (with <BinaryVolume> indices in `v`) to the specified (`r`, `g`, `b`, `a`) color.
 * @v the <BinaryVolume> indices of the voxels of interest
 * @pm the map from <BinaryVolume> indices to index of visualized point
 * @r the red value of the color for drawing the points in `v`
 * @g the green value of the color for drawing the points in `v`
 * @b the blue value of the color for drawing the points in `v`
 * @a the alpha value of the color for drawing the points in `v`
 *
 * Updates the colors of all the points specified in `v` by their <BinaryVolume> index to the specified (`r`, `g`, `b`, `a`) color.
 */
template <class T>
void setPointSetColor(const vector<T> &v, const map<voxelType, unsigned int> &pm,
		double r, double g, double b, double a){
	double c[4];
	c[0] = r;
	c[1] = g;
	c[2] = b;
	c[3] = a;
	map<unsigned int, vector<double> > uc;
	for(unsigned long i(0); i < v.size(); i++){
		if(pm.find((unsigned int)v[i]) == pm.end())
			continue;
		unsigned int pmIndex(pm.at((unsigned int)v[i]));
		for(unsigned int j(0); j < 4; j++)
			uc[pmIndex].push_back(c[j]);
	}
	wxGetApp().myFrame->m_canvas->updateCols(uc);
}

/* Sets all points to (`r`, `g`, `b`, `a`).
 * @r the red value of the color
 * @g the green value of the color
 * @b the blue value of the color
 * @a the opacity of the color
 *
 * Wrapper function to set all points to the color (`r`, `g`, `b`, `a`).
 */
void setAllPointsTo(double r, double g, double b, double a){
	wxGetApp().myFrame->m_canvas->updateCols(r, g, b, a);
}

/* Sets the point with <BinaryVolume> index `i` to the specified (`r`, `g`, `b`, `a`) color.
 * @i the <BinaryVolume> index of the voxel of interest
 * @pm the map from <BinaryVolume> indices to index of visualized point
 * @r the red value of the color for drawing the points in `v`
 * @g the green value of the color for drawing the points in `v`
 * @b the blue value of the color for drawing the points in `v`
 * @a the alpha value of the color for drawing the points in `v`
 *
 * Updates the color of the points specified by <BinaryVolume> index `i` to the specified (`r`, `g`, `b`, `a`) color.
 */
void setPointColor(voxelType i, const map<voxelType, unsigned int> &pm,
		double r, double g, double b, double a){
	if(pm.find(i) == pm.end())
		return;
	double c[4];
	c[0] = r;
	c[1] = g;
	c[2] = b;
	c[3] = a;
	map<unsigned int, vector<double> > uc;
	unsigned int pmIndex(pm.at(i));
	for(unsigned int j(0); j < 4; j++)
		uc[pmIndex].push_back(c[j]);
	wxGetApp().myFrame->m_canvas->updateCols(uc);
}

#endif
