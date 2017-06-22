#!/usr/bin/env python
# coding: UTF-8
#
## @package TextureMap
#
# Exercises texture mapping, by means of a
# textured cube, spinning around the X axis. 
#
# @author Paulo Cavalcanti
# @since 27/02/2017
#
# @see http://pyopengl.sourceforge.net/context/tutorials/nehe6.html
# @see https://open.gl/textures
#
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import time, sys, math
try:
    from PIL.Image import open
except ImportError as err:
    from Image import open

WIDTH = 640
HEIGHT = 480
ESCAPE = b'\033'

class TestContext( object ):
	"""NeHe 6 Demo"""

	initialPosition = (0,0,0) 

	## Set initial camera position, tutorial does the re-positioning
	#  OnInit is called by the Context class after initialization of the context has completed, 
	#  and before any rendering is attempted. Within this method, you'll generally perform 
	#  your global setup tasks.
	def __init__( self ):
		"""Load the image on initial load of the application"""

		self.imageID = self.loadImage ()

	## We are going to use the Python Imaging Library (PIL) for loading images, 
	#  something which is obviously not seen in the original tutorial.
	#  This method combines all of the functionality required to load the image with PIL, 
	#  convert it to a format compatible with PyOpenGL, generate the texture ID, 
	#  and store the image data under that texture ID.
	def loadImage( self, imageName = "images/velazquez_texture_256.jpg" ):
		"""Load an image file as a 2D texture using PIL"""

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

	def setCamera(self):
		"""Set camera parameters."""

		glMatrixMode(GL_PROJECTION)
		# Reset the Projection Matrix
		glLoadIdentity()

		# field of view, aspect ratio, near and far
		r = 1.0*math.sqrt(2.0) # ray of a sphere enclosing the object
		d = 6.0 # distance from camera to centre of the sphere (plane z = 0)
		fovy = (2.0*math.asin(r/d))*180/math.pi
		gluPerspective(fovy, float(WIDTH)/float(HEIGHT), 1.0, 10000.0)

		glMatrixMode (GL_MODELVIEW)
		glLoadIdentity()
		gluLookAt(3.5, 0, -6, 0, 0, 0, 0, 1, 0)

	def Render( self ):
		"""Render scene geometry"""

		self.setCamera()
		glRotated(time.time()%(8.0)/8 * -360, 1,0,0)

		self.setupTexture()
		self.drawCube()

	## The function called whenever a key is pressed.
	#  Note the use of Python tuples to pass in: (key, x, y)
	def keyPressed(self, *args):
		"""If escape is pressed, kill everything."""

		key = args [0]
		if key == ESCAPE:
			sys.exit ()

	## This method encapsulates the functions required to set up for textured rendering. 
	#  The original tutorial made these calls once for the entire program. 
	#  This organization makes more sense if you are likely to have multiple textures.
	def setupTexture( self ):
		"""Render-time texture environment setup"""

		# Configure the texture rendering parameters
		glEnable(GL_TEXTURE_2D)

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)

		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)
		# Re-select our texture, could use other generated textures if we had generated them earlier...
		glBindTexture(GL_TEXTURE_2D, imageID)

	## Drawing the cube has changed slightly, because we now need to specify the texture 
	#  coordinates for each vertex. This is all just taken from the original tutorial.
	def drawCube( self ):
		"""Draw a cube with texture coordinates"""

		# Clear Screen And Depth Buffer
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

		glBegin(GL_QUADS)
		glColor3f(1.0,0.0,0.0)
		glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, -1.0,  1.0);
		glTexCoord2f(1.0, 0.0); glVertex3f( 1.0, -1.0,  1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f( 1.0,  1.0,  1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f(-1.0,  1.0,  1.0);
		glColor3f(0.0,1.0,0.0)
		glTexCoord2f(1.0, 0.0); glVertex3f(-1.0, -1.0, -1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f(-1.0,  1.0, -1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f( 1.0,  1.0, -1.0);
		glTexCoord2f(0.0, 0.0); glVertex3f( 1.0, -1.0, -1.0);
		glColor3f(0.0,0.0,1.0)
		glTexCoord2f(0.0, 1.0); glVertex3f(-1.0,  1.0, -1.0);
		glTexCoord2f(0.0, 0.0); glVertex3f(-1.0,  1.0,  1.0);
		glTexCoord2f(1.0, 0.0); glVertex3f( 1.0,  1.0,  1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f( 1.0,  1.0, -1.0);
		glColor3f(1.0,1.0,0.0)
		glTexCoord2f(1.0, 1.0); glVertex3f(-1.0, -1.0, -1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f( 1.0, -1.0, -1.0);
		glTexCoord2f(0.0, 0.0); glVertex3f( 1.0, -1.0,  1.0);
		glTexCoord2f(1.0, 0.0); glVertex3f(-1.0, -1.0,  1.0);
		glColor3f(0.0,1.0,1.0)
		glTexCoord2f(1.0, 0.0); glVertex3f( 1.0, -1.0, -1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f( 1.0,  1.0, -1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f( 1.0,  1.0,  1.0);
		glTexCoord2f(0.0, 0.0); glVertex3f( 1.0, -1.0,  1.0);
		glColor3f(1.0,0.0,1.0)
		glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, -1.0, -1.0);
		glTexCoord2f(1.0, 0.0); glVertex3f(-1.0, -1.0,  1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f(-1.0,  1.0,  1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f(-1.0,  1.0, -1.0);
		glEnd()

		glDisable(GL_TEXTURE_2D)
		glutSwapBuffers()

	def OnIdle( self ):
		"""Request refresh of the context whenever idle"""

		self.Render()
		return 1

## Main program for testing.
def main(argv=None):
	global window

	if argv is None:
		argv = sys.argv

	# pass arguments to init
	glutInit(argv)

	# Select type of Display mode:   
	#  Double buffer 
	#  RGBA color
	#  Alpha components supported 
	#  Depth buffer
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
	
	# get a 640 x 480 window 
	glutInitWindowSize(WIDTH, HEIGHT)
	
	# the window starts at the upper left corner of the screen 
	glutInitWindowPosition(0, 0)
	
	# Okay, like the C version we retain the window id to use when closing, but for those of you new
	# to Python, remember this assignment would make the variable local and not global
	# if it weren't for the global declaration at the start of main.
	window = glutCreateWindow(b"TextureMap")

	# GLUT When mouse buttons are clicked in window
	#glutMouseFunc (Upon_Click)

	# GLUT When the mouse moves
	#glutMotionFunc (Upon_Drag)

	# We've told Glut the type of window we want, and we've told glut about
	# various functions that we want invoked (idle, resizing, keyboard events).
	# Glut has done the hard work of building up thw windows DC context and 
	# tying in a rendering context, so we are ready to start making immediate mode
	# GL calls.
	# Call to perform inital GL setup (the clear colors, enabling modes

	# This Will Clear The Background Color To Black
	glClearColor(0.0, 0.0, 0.0, 1.0)
	# Enables Clearing Of The Depth Buffer
	glClearDepth(1.0)
	# The Type Of Depth Test To Do
	glDepthFunc(GL_LEQUAL)
	# Enables Depth Testing
	glEnable(GL_DEPTH_TEST)
	# Select Flat Shading (Nice Definition Of Objects)
	glShadeModel (GL_SMOOTH)
	# Really Nice Perspective Calculations
	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)

	# mcolor will be applied to both ambient and diffuse components of the material.
	# This is done for convenience because in most cases Ambient and Diffuse properties
	# of a material should be set to the same color.
	mcolor = [ 1.0, 0.0, 0.0, 1.0 ]
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mcolor)
	glEnable (GL_LIGHT0)
	glEnable (GL_LIGHTING)
	# enable color tracking (specify material properties by merely calling the glColor)
	glEnable (GL_COLOR_MATERIAL)
	# set material properties which will be assigned by glColor
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE)
	glFrontFace(GL_CCW)

	# Create the GL context.
	glContext = TestContext()

	# Register the drawing function with glut, BUT in Python land, at least using PyOpenGL, we need to
	# set the function pointer and invoke a function to actually register the callback, otherwise it
	# would be very much like the C version of the code.	
	glutDisplayFunc(glContext.Render)
	
	# Register the function called when our window is resized.
	#glutReshapeFunc(glContext.reSizeScene)
	
	# Uncomment this line to get full screen.
	#glutFullScreen()

	# Register the function called when the keyboard is pressed.  
	glutKeyboardFunc(glContext.keyPressed)

	# When we are doing nothing, redraw the scene.
	glutIdleFunc(glContext.OnIdle)

	# Start Event Processing Engine	
	glutMainLoop()

## Print message to console, and kick off the main to get it rolling.
if __name__ == "__main__":
	print ("Hit ESC key to quit.")
	sys.exit(main())
