from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys
import copy
import random
from time import sleep
from math import cos, sin, degrees, acos
from plyfile import PlyData, PlyElement
from graph import *
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
plydata = PlyData.read('cube.ply')

pontos = plydata.elements[0].data
edges = plydata.elements[1].data
g = {}
selectedPolygonIndex = None
a = 0
polygons=[]

# Cria os poligonos e adiciona na lista
for edge in edges:
		for e in edge:
			points=[]
			for i in range(len(e)):
				points.append(Point(pontos[e[i]][0],pontos[e[i]][1],pontos[e[i]][2]))
			p = Polygon(points)
			p.points_indexes = e.tolist()
			polygons.append(p)

graph = Graph(g)
graph.generate_graph(g, edges, polygons)
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
	
	
	return True


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

	
	global g_isDragging, g_LastRot, g_Transform, g_ThisRot,polygons, line, selectedPolygonIndex

	g_isDragging = False
	if (button == GLUT_RIGHT_BUTTON and button_state == GLUT_UP):
		# Right button click
		# Obtem os pontos de clique do mouse
		xm,ym,zm = getMouse(cursor_x,cursor_y, 0)

		# Transforma os pontos de click do mouse
		p1 = dot(g_ThisRot,[xm,ym,-100])
		p2 = dot(g_ThisRot,[xm,ym,100])
		
		# Cria linha entre os dois pontos
		line = Line(Point(p1[0],p1[1],p1[2]),Point(p2[0],p2[1],p2[2]))

		### Pegar o poligono selecionado
		selectedPolygonIndex = False
		z = -999
		for index in range(len(polygons)):
			# Obtem dados sobre intercessao de um poligono na linha do click
			intersec = polygons[index].doesLineCrossPolygon(line)
			# Se a intercessao da linha for true
			# e o z do ponto interceptado for maior que o z anterior,
			# atualiza novo poligono e z mais na frente
			if intersec[0] and intersec[2] > z:
				selectedPolygonIndex = index
				z = intersec[2]

		# DEBUG: printar o poligono selecionado
		print "SELECTED FACE: " + str(selectedPolygonIndex)

		# g_LastRot = Matrix3fSetIdentity ();							# // Reset Rotation
		# g_ThisRot = Matrix3fSetIdentity ();							# // Reset Rotation
		# g_Transform = Matrix4fSetRotationFromMatrix3f (g_Transform, g_ThisRot);	# // Reset Rotation
	elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_UP):
		# Left button released
		g_LastRot = copy.copy (g_ThisRot);							# // Set Last Static Rotation To Last Dynamic One
	

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
	[1.0,  1.0,  0.0], [0.0,  0.502,  0.0], [0.0,  0.0,  1.0],
	[1.0, 0.0, 0.0], [1.0, 0.647, 0.0], [1.0, 1.0, 1.0],
	[1.0,  1.0,  0.0], [0.0,  0.502,  0.0], [0.0,  0.0,  1.0],
	[1.0, 0.0, 0.0], [1.0, 0.647, 0.0], [1.0, 1.0, 1.0],
	[1.0,  1.0,  0.0], [0.0,  0.502,  0.0], [0.0,  0.0,  1.0]];
	face = 0;

	for polygon in polygons:
		glBegin(GL_POLYGON);
		glColor3f(colors[face][0], colors[face][1], colors[face][2]);
		for point in polygon.points:
			glVertex3f(point.x,point.y,point.z)
		glEnd();
		face = face + 1;
	return

def rotateAndDraw(polygon, origen_polygon, rotate_point, rotate_axis, matrix=None):	
	global a
	b=0
	sense = rotate_axis.tripleProd(origen_polygon.normal,polygon.normal)
	n = len(polygons)
	dot_prod = origen_polygon.normal.dotProd(polygon.normal)
	# if dot_prod == 0.0:
	# 	angulo = 90
	# else:
	aux = dot_prod/origen_polygon.normal.len()*polygon.normal.len()
	angulo = math.degrees(math.acos(aux))
	angulo = round(angulo*(n-1))
	print angulo
	if a >= angulo:
		b = 0
	elif a >=0:
		b = -(0.5*sense)
		a += 0.5
	
	##### PREPARA MATRIZ DE TRANSFORMACAO
	## DEFINE PONTO DO EIXO DE ROTACAO TRANSLADADO
	##p = Point(1.0, 1.0, 1.0); ##Usando rotate_point
	## DEFINE QUAL SERA O EIXO DE ROTACAO
	##p_axis = Point(1.0, 0.0, 0.0); Usando rotate axis
	TR = translateAndRotate(b, rotate_point, rotate_axis)
	if matrix!=None:
		TR = dot(TR, matrix)

	polygon.matrix = TR
	##### TRANSFORMACOES APLICADAS EM CADA VERTICE
	for x in xrange(0, len(polygon.points)):
		aux_vertice_matrix = [polygon.points[x].x,polygon.points[x].y,polygon.points[x].z, 1]
		result_matrix = dot(TR, aux_vertice_matrix).tolist()[0] 
		polygon.points[x].x = result_matrix[0]
		polygon.points[x].y = result_matrix[1]
		polygon.points[x].z = result_matrix[2]

def rotateDede(root):	
	l = graph.breadth_first_search(root)
	
	for i in l['order']:
		parent = l['parent'][i]		
		aux_axis = polygons[parent].normal.crossProd(polygons[i].normal)
		aux2 = pontos[polygons[parent].edges[i][0]]
		aux_p = Point(aux2[0],aux2[1],aux2[2])
		if parent == root:
			rotateAndDraw(polygons[i], polygons[parent], aux_p, aux_axis)
		else:
			aux_axis = polygons[parent].normal.crossProd(polygons[i].normal)
			aux3 = polygons[parent].points_indexes.index(polygons[parent].edges[i][0])
			aux2 = polygons[parent].points[aux3]			
			aux_p = Point(aux2[0],aux2[1],aux2[2])
			rotateAndDraw(polygons[i], polygons[parent], aux_p, aux_axis, polygons[parent].matrix)
	DrawPolygon();

def Draw ():
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				# // Clear Screen And Depth Buffer
	glLoadIdentity();												# // Reset The Current Modelview Matrix
	glTranslatef(0.0,0.0,-10.0);									# // Move Left 1.5 Units And Into The Screen 6.0

	glPushMatrix();													# // NEW: Prepare Dynamic Transform
	glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform
	
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

	if selectedPolygonIndex == None:
		DrawPolygon()
	else:
		rotateDede(selectedPolygonIndex);

	glPopMatrix(); 												# // NEW: Unapply Dynamic Transform

	glLoadIdentity();												# // Reset The Current Modelview Matrix
	glTranslatef(0.0,0.0,-6.0);										# // Move Right 1.5 Units And Into The Screen 7.0

	glPushMatrix();													# // NEW: Prepare Dynamic Transform
	glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform

	glPopMatrix();	
	
	glFlush ();														# // Flush The GL Rendering Pipeline
	glutSwapBuffers()

	return