from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys
import copy
import random
from time import sleep
from math import cos, sin
from plyfile import PlyData, PlyElement
from graph import Graph
from matrix import *
from geometry import *

from ArcBall import * 				# ArcBallT and this tutorials set of points/vectors/matrix types

#PI2 = 2.0*3.1415926535			# 2 * PI (not squared!) 		// PI Squared

# *********************** Globals *********************** 
# Python 2.2 defines these directly
try:
	True
except NameError:
	True = 1==1
	False = 1==0

g_Transform = Matrix4fT ()
g_LastRot = Matrix3fT ()
g_ThisRot = Matrix3fT ()

g_ArcBall = ArcBallT (640, 480)
g_isDragging = False
g_quadratic = None
line = Line(Point(0,0,0.1),Point(0,0,0))
plydata = PlyData.read('dodecaedro.ply')
pontos = plydata.elements[0].data
edges = plydata.elements[1].data
g = {}
a = 0
polygons=[]

# Cria os poligonos e adiciona na lista
for edge in edges:
		for e in edge:
			points=[]
			for i in range(len(e)):
				points.append(Point(pontos[e[i]][0],pontos[e[i]][1],pontos[e[i]][2]))
			polygons.append(Polygon(points))


# A general OpenGL initialization function.  Sets all of the initial parameters. 
def Initialize (Width, Height):				# We call this right after our OpenGL window is created.
	global g_quadratic

	glClearColor(0.0, 0.0, 0.0, 1.0)					# This Will Clear The Background Color To Black
	glClearDepth(1.0)									# Enables Clearing Of The Depth Buffer
	glDepthFunc(GL_LEQUAL)								# The Type Of Depth Test To Do
	glEnable(GL_DEPTH_TEST)								# Enables Depth Testing
	glShadeModel (GL_FLAT);								# Select Flat Shading (Nice Definition Of Objects)
	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST) 	# Really Nice Perspective Calculations

	g_quadratic = gluNewQuadric();
	gluQuadricNormals(g_quadratic, GLU_SMOOTH);
	gluQuadricDrawStyle(g_quadratic, GLU_FILL); 
	# Why? this tutorial never maps any textures?! ? 
	# gluQuadricTexture(g_quadratic, GL_TRUE);			# // Create Texture Coords

	#glEnable (GL_LIGHT0)
	#glEnable (GL_LIGHTING)

	glEnable (GL_COLOR_MATERIAL)
	
	generate_graph()
	return True

def check_adjacent_face_has_vertices(faces, index, v1,v2):
    i = 0;
    for face in faces:
        if i != index:
            count = 0;
            for vertice in face:
                if vertice == v1 or vertice == v2:
                    count = count + 1
            if count == 2:
                return i;
        i = i+1;

def generate_graph():
	i = 0;
	faces = [];
	for edge in edges :
	    faces.append(edges[i][0].tolist())
	    g[i] = [] 
	    i = i+1
	
	graph = Graph(g)
	i=0;
	for face in faces:
	    for vertice in range(len(face)-1):
	        viz = check_adjacent_face_has_vertices(faces,i,face[vertice],face[vertice+1])
	        graph.add_edge({i,viz})
	    viz = check_adjacent_face_has_vertices(faces,i,face[0],face[vertice+1])
	    graph.add_edge({i,viz})   
	    i=i+1
	#print(g)    
	return

def getMouse(cursor_x, cursor_y, z):
	modelView = glGetDoublev( GL_MODELVIEW_MATRIX );
	projection = glGetDoublev( GL_PROJECTION_MATRIX );
	viewport = glGetIntegerv( GL_VIEWPORT );
	cursor_x = float (cursor_x);
	cursor_y = float (viewport[3]) - float (cursor_y);
	cursor_z = glReadPixels(cursor_x, cursor_y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT);
	posX, posY, posZ = gluUnProject(cursor_x, cursor_y, cursor_z, modelView, projection, viewport);
	return posX, posY, posZ;


def Upon_Drag (cursor_x, cursor_y):
	""" Mouse cursor is moving
		Glut calls this function (when mouse button is down)
		and pases the mouse cursor postion in window coords as the mouse moves.
	"""
	
	global g_isDragging, g_LastRot, g_Transform, g_ThisRot

	if (g_isDragging):
		mouse_pt = Point2fT (cursor_x, cursor_y)
		ThisQuat = g_ArcBall.drag (mouse_pt)						# // Update End Vector And Get Rotation As Quaternion
		g_ThisRot = Matrix3fSetRotationFromQuat4f (ThisQuat)		# // Convert Quaternion Into Matrix3fT
		# Use correct Linear Algebra matrix multiplication C = A * B
		g_ThisRot = Matrix3fMulMatrix3f (g_LastRot, g_ThisRot)		# // Accumulate Last Rotation Into This One
		g_Transform = Matrix4fSetRotationFromMatrix3f (g_Transform, g_ThisRot)	# // Set Our Final Transform's Rotation From This One
	return

def Upon_Click (button, button_state, cursor_x, cursor_y):
	""" Mouse button clicked.
		Glut calls this function when a mouse button is
		clicked or released.
	"""

	
	global g_isDragging, g_LastRot, g_Transform, g_ThisRot,polygons, line

	g_isDragging = False
	if (button == GLUT_RIGHT_BUTTON and button_state == GLUT_UP):
		# Right button click
		g_LastRot = Matrix3fSetIdentity ();							# // Reset Rotation
		g_ThisRot = Matrix3fSetIdentity ();							# // Reset Rotation
		g_Transform = Matrix4fSetRotationFromMatrix3f (g_Transform, g_ThisRot);	# // Reset Rotation
	elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_UP):
		# Left button released
		g_LastRot = copy.copy (g_ThisRot);							# // Set Last Static Rotation To Last Dynamic One


		# Obtem os pontos de clique do mouse
		xm,ym,zm = getMouse(cursor_x,cursor_y, 0)

		# Transforma os pontos de click do mouse
		p1 = dot(g_ThisRot,[xm,ym,-100])
		p2 = dot(g_ThisRot,[xm,ym,100])
		
		# Cria linha entre os dois pontos
		line = Line(Point(p1[0],p1[1],p1[2]),Point(p2[0],p2[1],p2[2]))

		### Pegar o poligono selecionado
		selectedPolygon = False
		z = -999
		for polygon in polygons:
			# Obtem dados sobre intercessao de um poligono na linha do click
			intersec = polygon.doesLineCrossPolygon(line)
			# Se a intercessao da linha for true
			# e o z do ponto interceptado for maior que o z anterior,
			# atualiza novo poligono e z mais na frente
			if intersec[0] and intersec[1][2] > z:
				selectedPolygon = polygon
				z = intersec[1][2]

		# DEBUG: printar o poligono selecionado
		print selectedPolygon

	elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_DOWN):	
		# Left button clicked down
		g_LastRot = copy.copy (g_ThisRot);							# // Set Last Static Rotation To Last Dynamic One
		g_isDragging = True											# // Prepare For Dragging
		mouse_pt = Point2fT (cursor_x, cursor_y)
		g_ArcBall.click (mouse_pt);								# // Update Start Vector And Prepare For Dragging

	return


def DrawPolygon():
	colors = [[1.0, 0.0, 0.0], [1.0, 0.647, 0.0], [1.0, 1.0, 1.0],
	[1.0,  1.0,  0.0], [0.0,  0.502,  0.0], [0.0,  0.0,  1.0],
	[1.0, 0.0, 0.0], [1.0, 0.647, 0.0], [1.0, 1.0, 1.0],
	[1.0,  1.0,  0.0], [0.0,  0.502,  0.0], [0.0,  0.0,  1.0]];
	face = 0;
	# for edge in edges:
		# glBegin(GL_POLYGON);
		# glColor3f(colors[face][0], colors[face][1], colors[face][2]);
	# 	for e in edge:
	# 		for i in range(len(e)):
	# 			glVertex3f(pontos[e[i]][0], pontos[e[i]][1],  pontos[e[i]][2]);

	# 	glEnd();
	for polygon in polygons:
		glBegin(GL_POLYGON);
		glColor3f(colors[face][0], colors[face][1], colors[face][2]);
		for point in polygon.points:
			glVertex3f(point.x,point.y,point.z)
		glEnd();
		face = face + 1;
	return


def Cube():	
	#FRENTE	
	glBegin(GL_POLYGON); 
	glColor3f(   1.0,  0.0, 0.0 );
	glVertex3f(  0.5, -0.5, 0.5 );     
	glVertex3f(  0.5,  0.5, 0.5 );      
	glVertex3f( -0.5,  0.5, 0.5 );      
	glVertex3f( -0.5, -0.5, 0.5 );       
	glEnd();

	#TRASEIRA
	glBegin(GL_POLYGON);
	glColor3f(   1.0,  0.647, 0.0 );
	glVertex3f(  0.5, -0.5, -0.5 );
	glVertex3f(  0.5,  0.5, -0.5 );
	glVertex3f( -0.5,  0.5, -0.5 );
	glVertex3f( -0.5, -0.5, -0.5 );
	glEnd();
	 
	#DIREITA
	glBegin(GL_POLYGON);
	glColor3f(  1.0,  1.0,  1.0 );
	glVertex3f( 0.5, -0.5, -0.5 );
	glVertex3f( 0.5,  0.5, -0.5 );
	glVertex3f( 0.5,  0.5,  0.5 );
	glVertex3f( 0.5, -0.5,  0.5 );
	glEnd();
	 
	#ESQUERDA
	glBegin(GL_POLYGON);
	glColor3f(   1.0,  1.0,  0.0 );
	glVertex3f( -0.5, -0.5,  0.5 );
	glVertex3f( -0.5,  0.5,  0.5 );
	glVertex3f( -0.5,  0.5, -0.5 );
	glVertex3f( -0.5, -0.5, -0.5 );
	glEnd();
	 
	#TOPO
	glBegin(GL_POLYGON);
	glColor3f(   0.0,  0.502,  0.0 );
	glVertex3f(  0.5,  0.5,  0.5 );
	glVertex3f(  0.5,  0.5, -0.5 );
	glVertex3f( -0.5,  0.5, -0.5 );
	glVertex3f( -0.5,  0.5,  0.5 );
	glEnd();
	 
	#BASE
	glBegin(GL_POLYGON);
	glColor3f(   0.0,  0.0,  1.0 );
	glVertex3f(  0.5, -0.5, -0.5 );
	glVertex3f(  0.5, -0.5,  0.5 );
	glVertex3f( -0.5, -0.5,  0.5 );
	glVertex3f( -0.5, -0.5, -0.5 );
	glEnd();														# // Done Torus
	return

def josh():
	## JOSH MODE ON

	##### MATRIZ DO POLIGONO
	polygons = [[ 0.5, -0.5, 0.5, 1],
			[ 0.5,  0.5, 0.5, 1],
			[-0.5,  0.5, 0.5, 1],
			[-0.5, -0.5, 0.5, 1]];
	

	# ##### POLIGONO VERDE INICIAL
	glBegin(GL_POLYGON);
	glColor3f( 0.0, 1.0, 0.0 ); # verde
	for line in polygons:
		glVertex3f(  line[1-1], line[2-1], line[3-1] );     
	glEnd(); 


	##### PREPARA MATRIZ DE TRANSFORMACAO
	## DEFINE PONTO DO EIXO DE ROTACAO TRANSLADADO
	p = Point(0.0, -0.5, 0.5);
	## DEFINE QUAL SERA O EIXO DE ROTACAO
	p_axis = Point(1.0, 0.0, 0.0);
	TR = translateAndRotate(90, p, p_axis)


	##### TRANSFORMACOES APLICADAS EM CADA VERTICE
	polygons[0] = dot(TR, polygons[0]).tolist()[0]
	polygons[1] = dot(TR, polygons[1]).tolist()[0]
	polygons[2] = dot(TR, polygons[2]).tolist()[0]
	polygons[3] = dot(TR, polygons[3]).tolist()[0]

	
	##### POLIGONO VERMELHO TRANSFORMADO
	glBegin(GL_POLYGON);
	glColor3f( 1.0, 0.0, 0.0 ); # vermelho
	for line in polygons:
		glVertex3f(  line[1-1], line[2-1], line[3-1] );     
	glEnd(); 
	## JOSH MODE OFF
	return

def rotateAndDraw(polygon):
	basic_axis = {'x': Point(1.0,0.0,0.0), 'y': Point(0.0,1.0,0.0), 'z': Point(0.0,0.0,1.0)}
	global a
	b=0
	if a > 90:
		b = 0
		print(polygon)
		#a -= 2
	elif a >=0:
		b = 2
		a += 2

	
	##### PREPARA MATRIZ DE TRANSFORMACAO
	## DEFINE PONTO DO EIXO DE ROTACAO TRANSLADADO
	p = Point(1.0, 1.0, 1.0); ##Usando rotate_point
	## DEFINE QUAL SERA O EIXO DE ROTACAO
	##p_axis = Point(1.0, 0.0, 0.0); Usando rotate axis
	TR = translateAndRotate(b, p, basic_axis['x'])
	#polygons[0] = Polygon([Point(1.0,1.0,1.0),Point(2.0,2.0,2.0),Point(3.0,3.0,3.0)])

	##### TRANSFORMACOES APLICADAS EM CADA VERTICE
	#print(polygon.points)
	
	# aux_vertice_matrix = [polygon.points[0][0],polygon.points[0][1],polygon.points[0][2], 1]
	# result_matrix = dot(TR, aux).tolist()[0] 
	# print(result_matrix);
	# exit();
	for x in xrange(0, len(polygon.points)):
		aux_vertice_matrix = [polygon.points[x][0],polygon.points[x][1],polygon.points[x][2], 1]
		result_matrix = dot(TR, aux_vertice_matrix).tolist()[0] 
		polygon.points[x].x = result_matrix[0]
		polygon.points[x].y = result_matrix[1]
		polygon.points[x].z = result_matrix[2]

	# print(polygon.points)
	# exit();
	DrawPolygon();
	#exit();

def rotateDede():
	basic_axis = {'x': Point(1.0,0.0,0.0), 'y': Point(0.0,1.0,0.0), 'z': Point(0.0,0.0,1.0)}

	for polygon in polygons:
		glBegin(GL_POLYGON);
		glColor3f(colors[face][0], colors[face][1], colors[face][2]);
		for point in polygon.points:
			glVertex3f(point.x,point.y,point.z)
		glEnd();
		face = face + 1;
	return

def Draw ():
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				# // Clear Screen And Depth Buffer
	glLoadIdentity();												# // Reset The Current Modelview Matrix
	glTranslatef(0.0,0.0,-6.0);									# // Move Left 1.5 Units And Into The Screen 6.0

	glPushMatrix();													# // NEW: Prepare Dynamic Transform
	glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform
	# glColor3f(0.0,0.0,0.0);
	# DrawPolygon();
	#Cube();

	### Eixo global x vermelho
	glBegin(GL_LINES)
	glColor3f(1.0,0.0,0.0)
	glVertex2f(0.0, 0.0)
	glVertex2f(2.0, 0.0)
	glEnd()

	### Eixo global y amarelo
	glBegin(GL_LINES)
	glColor3f(1.0,1.0,0.0)
	glVertex2f(0.0, 0.0)
	glVertex2f(0.0, 2.0)
	glEnd()

	### Eixo global z magenta
	glBegin(GL_LINES)
	glColor3f(1.0,0.0,1.0)
	glVertex3f(0.0, 0.0, 0.0)
	glVertex3f(0.0, 0.0, 2.0)
	glEnd()

	### Linha de pick (line global)
	glBegin(GL_LINES)
	glColor3f(1.0,1.0,1.0)
	glVertex3f(line.p1[0], line.p1[1], line.p1[2])
	glVertex3f(line.p2[0], line.p2[1], line.p2[2])
	glEnd()



	### Normal do primeiro poligono
	polygon_test = polygons[0]
	glBegin(GL_LINES)
	glColor3f(0.0,1.0,1.0)
	glVertex3f(0.0, 0.0, 0.0)
	glVertex3f(-polygon_test.compNormal()[0], -polygon_test.compNormal()[1], -polygon_test.compNormal()[2])
	glEnd()



	#r = [0.5 - 0.5, 0.5 + 0.5, 0.5 - 0.5]
	# glRotatef(90,1.0,0.0,0.0);
	# glRotatef(90,0.5,-0.5,0.5);
	# glTranslatef(1,0,0)
	# glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
	
	DrawPolygon()
	# rotateAndDraw(polygons[4]);

	glPopMatrix(); 												# // NEW: Unapply Dynamic Transform

	glLoadIdentity();												# // Reset The Current Modelview Matrix
	glTranslatef(0.0,0.0,-6.0);										# // Move Right 1.5 Units And Into The Screen 7.0

	glPushMatrix();													# // NEW: Prepare Dynamic Transform
	glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform

	#glBegin(GL_POLYGON); 
	#glColor3f(   1.0,  0.0, 0.0 );
	#glVertex3f(  0.5, -0.5, 0.5 );     
	#glVertex3f(  0.5,  0.5, 0.5 );      
	#glVertex3f( -0.5,  0.5, 0.5 );      
	#glVertex3f( -0.5, -0.5, 0.5 );       
	#glEnd();
	#glColor3f(1.0,0.75,0.75);
	#gluSphere(g_quadratic,1.3,20,20);
	glPopMatrix();	
	
	glFlush ();														# // Flush The GL Rendering Pipeline
	glutSwapBuffers()

	return