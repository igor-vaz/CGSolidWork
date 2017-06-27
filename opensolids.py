from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys
import copy
import random as rd
from time import sleep
from math import cos, sin, degrees, acos
from plyfile import PlyData, PlyElement
from graph import *
from matrix import *
from geometry import *
import copy
from numpy.linalg import inv

from ArcBall import * 				# ArcBallT and this tutorials set of points/vectors/matrix types
try:
    from PIL.Image import open
except ImportError as err:
    from Image import open

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


if len(sys.argv) == 1:
	print "Use: python "+sys.argv[0]+" [modelo] [velocidade = 2] [exibir_eixos = False]"
	print "Modelos disponiveis:"
	print " - tetraedro"
	print " - octaedro"
	print " - hexaedro"
	print " - icosaedro"
	print " - dodecaedro"
	print " bonus:"
	print " - dart"
	print " - cow"
	print " - porsche"
	exit(1)

vertexs = None
edges = None
polygons=[]
colors=[]
graph = None
selectedFace = False
animateProgress = None
zoom = -8
imageID = None
isSolidOpen = False

# A general OpenGL initialization function.  Sets all of the initial parameters. 
def Initialize (Width, Height):				# We call this right after our OpenGL window is created.
	global g_quadratic, vertexs, edges, polygons, graph, animateProgress, colors, imageID

	glClearColor(0.0, 0.0, 0.0, 1.0)					# This Will Clear The Background Color To Black
	glClearDepth(1.0)									# Enables Clearing Of The Depth Buffer
	glDepthFunc(GL_LEQUAL)								# The Type Of Depth Test To Do
	glEnable(GL_DEPTH_TEST)								# Enables Depth Testing
	glShadeModel (GL_FLAT);								# Select Flat Shading (Nice Definition Of Objects)
	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST) 	# Really Nice Perspective Calculations

	g_quadratic = gluNewQuadric();
	gluQuadricNormals(g_quadratic, GLU_SMOOTH);
	gluQuadricDrawStyle(g_quadratic, GLU_FILL); 
	glEnable (GL_COLOR_MATERIAL)
	
	# Obtem modelo e separa os vertices e bordas
	plydata = PlyData.read('modelos/'+sys.argv[1]+'.ply')
	vertexs = plydata.elements[0].data
	edges = plydata.elements[1].data
	
	# Cria os poligonos e adiciona na lista de poligonos
	for edge in edges:
		edge = edge[0]
		points=[]
		for i in range(len(edge)):
			points.append(Point(vertexs[edge[i]][0],vertexs[edge[i]][1],vertexs[edge[i]][2]))
		polygon = Polygon(points)
		polygon.points_indexes = edge.tolist()
		# Guarda a normal inicial antes de sofrer as transformacoes
		polygon.original_normal = polygon.normal
		polygons.append(polygon)
	
	for i in range(len(polygons)):
		colors.append([rd.random(),rd.random(),rd.random()])

	# Cria o grafo e gera baseado na lista de poligonos e bordas
	graph = Graph({})
	graph.generate_graph({}, edges, polygons)
	
	animateProgress = [0] * len(polygons)

	imageID = loadImage()

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

	global g_isDragging, g_LastRot, g_Transform, g_ThisRot,polygons, selectedFace, zoom



	if button == 1 and button_state == GLUT_UP:
		flatMapSize()



	# Mouse Scroll
	if button == 3:
		zoom += float(button_state)/4
	if button == 4:
		zoom -= float(button_state)/4

	# Mouse click
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
		selectedFace = False
		z = -999
		for index in range(len(polygons)):
			# Obtem dados sobre intercessao de um poligono na linha do click
			intersec = polygons[index].doesLineCrossPolygon(line)
			# Se a intercessao da linha for true
			# e o z do ponto interceptado for maior que o z anterior,
			# atualiza novo poligono e z mais na frente
			if intersec[0] and intersec[2] > z:
				selectedFace = index
				z = intersec[2]

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
	global colors

	face = 0;
	for polygon in polygons:
		glBegin(GL_POLYGON);
		glColor3f(colors[face][0], colors[face][1], colors[face][2]);
		p = 0
		for point in polygon.points:
			if len(polygon.texture_coords) > 0:
				glTexCoord2f(polygon.texture_coords[p][0], polygon.texture_coords[p][1]);
			glVertex3f(point.x,point.y,point.z)
			p += 1
		glEnd();
		face = face + 1;
		if face >= len(colors):
			face = 0
	return

def rotateFace(polygon, polygon_origin, index, rotate_point, rotate_axis, matrixTransParent = None):	
	global isSolidOpen, animateProgress,selectedFace
	sense = rotate_axis.tripleProd(polygon_origin.original_normal,polygon.original_normal)

	# Produto vetorial entre as normas
	dot_prod = polygon_origin.original_normal.dotProd(polygon.original_normal)
	aux = dot_prod/polygon_origin.original_normal.len()*polygon.original_normal.len()
	# Define angulo de abertura
	angulo = math.degrees(math.acos(aux))
	
	### DEFINE VELOCIDADE DE ABERTURA SE PASSADO SEGUNDO ARGUMENTO
	speed = 2
	if len(sys.argv) >= 3:
		speed = float(sys.argv[2])

	### PREPARA SETORES DE b PARA ROTACIONAR FRACIONADO ATE O ANGULO FINAL COMO SE FOSSE ANIMADO
	b=0
	# if animateProgress[index] >= angulo:
	# 	b = 0
	# 	isSolidOpen = True
	# elif animateProgress[index] >=0:
	# 	b = -(speed)
	# 	print(b)
	# 	animateProgress[index] += speed
	
	#animacao para abrir o solido
	if animateProgress[index] < angulo and not isSolidOpen:
		b = -(speed)
		animateProgress[index] += speed
	#animacao para fechar o solido(EM TEORIA)
	elif animateProgress[index] > 0 and isSolidOpen:
		b = -speed
		animateProgress[index] -= speed
	#Solido fechado e face selecionada apagada
	# elif animateProgress[index] < speed and isSolidOpen:
	# 	isSolidOpen = False
	# 	selectedFace = False
	#Solido aberto
	elif animateProgress[index] >= angulo:
		b = 0
		isSolidOpen = True


	print(animateProgress[index])
	print(isSolidOpen)
	print("b = "+str(b))
	### PREPARA MATRIZ DE TRANSFORMACAO
	matrizTrans = translateAndRotate(b, rotate_point, rotate_axis)
	if matrixTransParent is not None and not isSolidOpen:
		matrizTrans = dot(matrizTrans, matrixTransParent)
	elif matrixTransParent is not None and isSolidOpen:
		matrizTrans = dot(matrizTrans, inv(matrixTransParent))

	polygon.matrix = matrizTrans
	### TRANSFORMACOES APLICADAS EM CADA VERTICE
	for x in xrange(0, len(polygon.points)):
		vertice_matrix = [polygon.points[x].x,polygon.points[x].y,polygon.points[x].z, 1]
		result_matrix = dot(matrizTrans, vertice_matrix).tolist()[0] 
		polygon.points[x].x = result_matrix[0]
		polygon.points[x].y = result_matrix[1]
		polygon.points[x].z = result_matrix[2]
	polygon.normal = polygon.compNormal().normalize()

def openFrom(root):
	tree = graph.breadth_first_search(root)
	
	for node in tree['order']:
		parent = tree['parent'][node]	
		p_axis = polygons[parent].normal.crossProd(polygons[node].normal)
		if parent == root:
			vertex = vertexs[polygons[parent].edges[node][0]]
			p_ref = Point(vertex[0],vertex[1],vertex[2])
			rotateFace(polygons[node], polygons[parent], node, p_ref, p_axis)
		else:
			aux = polygons[parent].points_indexes.index(polygons[parent].edges[node][0])
			vertex = polygons[parent].points[aux]			
			p_ref = Point(vertex[0],vertex[1],vertex[2])
			rotateFace(polygons[node], polygons[parent],node, p_ref, p_axis, polygons[parent].matrix)

def closeFrom(root):
	tree = graph.breadth_first_search(root)
	for node in tree['order']:
		parent = tree['parent'][node]		
		p_axis = polygons[parent].normal.crossProd(polygons[node].normal)
		if parent == root:
			vertex = vertexs[polygons[parent].edges[node][0]]
			p_ref = Point(vertex[0],vertex[1],vertex[2])
			rotateFace(polygons[node], polygons[parent], node, p_ref, p_axis)
		else:
			aux = polygons[parent].points_indexes.index(polygons[parent].edges[node][0])
			vertex = polygons[parent].points[aux]			
			p_ref = Point(vertex[0],vertex[1],vertex[2])
			rotateFace(polygons[node], polygons[parent],node, p_ref, p_axis, polygons[parent].matrix)

def closeFace(polygon, polygon_origin, index, rotate_point, rotate_axis, matrixTransParent = None):	
	global isSolidOpen,animateProgress,selectedFace
	sense = rotate_axis.tripleProd(polygon_origin.original_normal,polygon.original_normal)
	
	# Produto vetorial entre as normas
	dot_prod = polygon_origin.original_normal.dotProd(polygon.original_normal)
	aux = dot_prod/polygon_origin.original_normal.len()*polygon.original_normal.len()
	# Define angulo de abertura
	angulo = math.degrees(math.acos(aux))
	
	### DEFINE VELOCIDADE DE ABERTURA SE PASSADO SEGUNDO ARGUMENTO
	speed = 2
	if len(sys.argv) >= 3:
		speed = float(sys.argv[2])

	### PREPARA SETORES DE b PARA ROTACIONAR FRACIONADO ATE O ANGULO FINAL COMO SE FOSSE ANIMADO
	b=0
	if animateProgress[index] <= 0:
		b = 0
		isSolidOpen = False
		selectedFace = False
	elif animateProgress[index] >= angulo:
		b = -(speed)
		animateProgress[index] -= speed
	
	### PREPARA MATRIZ DE TRANSFORMACAO
	matrizTrans = translateAndRotate(b, rotate_point, rotate_axis)
	if matrixTransParent is not None:
		matrizTrans = dot(matrizTrans, matrixTransParent)

	polygon.matrix = matrizTrans
	### TRANSFORMACOES APLICADAS EM CADA VERTICE
	for x in xrange(0, len(polygon.points)):
		vertice_matrix = [polygon.points[x].x,polygon.points[x].y,polygon.points[x].z, 1]
		result_matrix = dot(matrizTrans, vertice_matrix).tolist()[0] 
		polygon.points[x].x = result_matrix[0]
		polygon.points[x].y = result_matrix[1]
		polygon.points[x].z = result_matrix[2]
	polygon.normal = polygon.compNormal().normalize()

def Draw ():
	global zoom, isSolidOpen,selectedFace

	if len(polygons[0].texture_coords) > 0:
		setupTexture()

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				# // Clear Screen And Depth Buffer
	glLoadIdentity();		

	### DEFINE ZOOM dinamico
	glTranslatef(0.0,0.0,zoom);

	glPushMatrix();													# // NEW: Prepare Dynamic Transform
	glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform

	### DEFINE EIXOS PRINCIPAIS SE PASSADO TERCEIRO ARGUMENTO
	if len(sys.argv) >= 4 and sys.argv[3] == "True":
		### Eixo global x
		glBegin(GL_LINES)
		glColor3f(1.0,0.0,0.0)
		glVertex2f(0.0, 0.0)
		glVertex2f(2.0, 0.0)
		glEnd()
		### Eixo global y
		glBegin(GL_LINES)
		glColor3f(1.0,1.0,0.0)
		glVertex2f(0.0, 0.0)
		glVertex2f(0.0, 2.0)
		glEnd()
		### Eixo global z
		glBegin(GL_LINES)
		glColor3f(1.0,0.0,1.0)
		glVertex3f(0.0, 0.0, 0.0)
		glVertex3f(0.0, 0.0, 2.0)
		glEnd()

	### ABRE SOLIDO SE TIVER POLIGONO SELECIONADO
	if selectedFace is not False:
		openFrom(selectedFace);

	### DESENHA POLIGONO
	DrawPolygon();

	if isSolidOpen:
		flatMapSize() 
		if selectedFace is not False:
			closeFrom(selectedFace)
		# isSolidOpen = False
		# selectedFace = False

	glDisable(GL_TEXTURE_2D)

	glPopMatrix();													# // NEW: Unapply Dynamic Transform
	glFlush ();														# // Flush The GL Rendering Pipeline
	glutSwapBuffers()

	return

## We are going to use the Python Imaging Library (PIL) for loading images, 
#  something which is obviously not seen in the original tutorial.
#  This method combines all of the functionality required to load the image with PIL, 
#  convert it to a format compatible with PyOpenGL, generate the texture ID, 
#  and store the image data under that texture ID.
def loadImage():
	"""Load an image file as a 2D texture using PIL"""
	imageName = "images/velazquez_texture_256.jpg" 
	# PIL defines an "open" method which is Image specific!
	im = open(imageName)
	try:
		ix, iy, image = im.size[0], im.size[1], im.tobytes("raw", "RGBA", 0, -1)
	except (SystemError, ValueError):
		ix, iy, image = im.size[0], im.size[1], im.tobytes("raw", "RGBX", 0, -1)
	except AttributeError:
		ix, iy, image = im.size[0], im.size[1], im.tostring("raw", "RGBX", 0, -1)

	# Generate a texture ID
	ID = glGenTextures(1)

	# Make our new texture ID the current 2D texture
	glBindTexture(GL_TEXTURE_2D, ID)
	glPixelStorei(GL_UNPACK_ALIGNMENT,1)

	# Copy the texture data into the current texture ID
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ix, iy, 0, GL_RGBA, GL_UNSIGNED_BYTE, image)

	# Note that only the ID is returned, no reference to the image object or the 
	# string data is stored in user space. 
	# The data is only present within the GL after this call exits.
	return ID

def setupTexture():
	global imageID
	"""Render-time texture environment setup"""

	# Configure the texture rendering parameters
	glEnable(GL_TEXTURE_2D)

	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)
	# Re-select our texture, could use other generated textures if we had generated them earlier...
	glBindTexture(GL_TEXTURE_2D, imageID)

def flatMapSize():
	global polygons

	polygons2 = copy.deepcopy(polygons)


	maxX, minX, maxY, minY = 0,0,0,0;

	for polygon in polygons2:
		selectedPolygon = polygon
		zAxis = Point(0,0,1)
		# Produto vetorial entre as normas
		dot_prod = selectedPolygon.normal.dotProd(zAxis)
		aux = dot_prod/selectedPolygon.normal.len()*zAxis.len()
		# Define angulo de abertura
		angulo = math.degrees(math.acos(aux))

		if angulo != 0 and angulo != 180:
			rotate_axis = selectedPolygon.normal.crossProd(zAxis)

			matrizTrans = translateAndRotate(angulo, Point(0,0.0,0), rotate_axis)

		for i in xrange(0, len(polygon.points)):
			if angulo != 0 and angulo != 180:
				vertice_matrix = [polygon.points[i].x,polygon.points[i].y,polygon.points[i].z, 1]
				result_matrix = dot(matrizTrans, vertice_matrix).tolist()[0] 
				polygon.points[i].x = result_matrix[0]
				polygon.points[i].y = result_matrix[1]
				polygon.points[i].z = result_matrix[2]

			point = polygon.points[i]

			if maxX < point.x: 
				maxX = point.x

			if minX > point.x: 
				minX = point.x 

			if maxY < point.y: 
				maxY = point.y

			if minY > point.y: 
				minY = point.y 

	largura = maxX - minX
	altura = maxY - minY

	for j in xrange(0, len(polygons2)):
		polygon = polygons2[j]
		for i in xrange(0, len(polygon.points)):
			point = polygon.points[i]

			mapX = (point.x - minX) / largura
			mapY = (point.y - minY) / altura

			polygons[j].texture_coords.append([mapX, mapY])

	return