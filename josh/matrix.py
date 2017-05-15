#!/usr/bin/env python
# coding: UTF-8
#
## @package matrix
#
# OpenGL is column major and numpy row major.
# We use row major here, but could add order='F' in reshape,
# to transform to column major, I guess...
#
# @author Paulo Cavalcanti
# @since 13/02/2017
# @see http://3dengine.org/Rotate_arb
# @see http://www.python-course.eu/matrix_arithmetic.php
#
import sys, math
import numpy as np 
from math import cos,sin

## Returns a translation matrix.
#
#  @param dx x translation.
#  @param dy y translation.
#  @param dz z translation.
#  @return translation matrix.
#
def translate(dx,dy,dz):
	M = identity()

	M[0,3] = dx
	M[1,3] = dy
	M[2,3] = dz

	return M

## Returns a scale matrix.
#
#  @param sx x scale.
#  @param sy y scale.
#  @param sz z scale.
#  @return scale matrix.
#
def scale(sx,sy,sz):
	M = identity()

	M[0,0] = sx
	M[1,1] = sy
	M[2,2] = sz

	return M

## Returns a rotation matrix.
#
#  @param ang rotation angle in degrees.
#  @param x rotation axis (x vector component).
#  @param y rotation axis (y vector component).
#  @param z rotation axis (z vector component).
#  @return rotation matrix.
#
def rotate(ang, x,y,z):
    ang = np.deg2rad(ang)
    c=cos(ang)
    s=sin(ang)
    t=1-c

    len = math.sqrt(x * x + y * y + z * z)
    len = 1 / len
    x *= len
    y *= len
    z *= len

    M = np.matrix([
        [t*x*x+c,    t*x*y-s*z,  t*x*z+s*y,  0],
        [t*x*y+s*z,  t*y*y+c,    t*y*z-s*x,  0],
        [t*x*z-s*y,  t*y*z+s*x,  t*z*z+c,    0],
        [0,  0,  0,  1],
    ])

    return M

## Returns an identity matrix.
#
#  Same as:
#    glPushMatrix()
#    glLoadIdentity()
#    c = glGetDoublev ( GL_MODELVIEW_MATRIX )
#    glPopMatrix()
#    return c
#
#  @return identity matrix.
#
def identity():
    M = np.matrix(np.identity(4))

    return M

## Matrix multiplication.
#
#  Same as:
#    glPushMatrix()
#    glLoadMatrixf(a)
#    glMultMatrixf(b)
#    c = glGetDoublev ( GL_MODELVIEW_MATRIX )
#    glPopMatrix()
#    return c
#
#  @param a first matrix.
#  @param b second matrix.
#  @return a x b.
#
def dot(a,b):
	return np.dot(a,b)

## Rotate around an axis, passing through a given point.
#
#   Same as:
#     glPushMatrix()
#     glLoadIdentity()
#     glTranslate(p.x,p.y,p.z)
#     glRotate(ang, axix.x, axis.y, axis.z)
#     glTranslate(-p.x,-p.y,-p.z)
#     T = glGetDoublev ( GL_MODELVIEW_MATRIX )
#     glPopMatrix()
#     return T 
#
#   @param ang rotation angle.
#   @param p point the axix passes through.
#   @param axis rotation axis.
#   @return rotation matrix: T.R.-T
#
def translateAndRotate(ang, p, axis):
	T = translate(p.x,p.y,p.z) * \
        rotate(ang, axis.x, axis.y, axis.z) * \
        translate(-p.x,-p.y,-p.z)
	return T

## Apply a given transformation t, using p as the fixed point.
def translateAndTransform(t, p):
	T = translate(p.x,p.y,p.z) * t * translate(-p.x,-p.y,-p.z)
	return T

## Return a rotation matrix, given three angles in the order: ZYX
#
#  @see http://www.chrobotics.com/library/understanding-euler-angles
#  @param angles a list with angle x, y and z.
#  @return a matrix Z * Y * X
#
def rotateZYX(angles):
	m1 = rotate(angles[2], 0.0, 0.0, 1.0)
	m2 = rotate(angles[1], 0.0, 1.0, 0.0)
	m3 = rotate(angles[0], 1.0, 0.0, 0.0)
	# m1 x m2 x m3
	m = dot(dot(m1,m2),m3)
	return np.asarray(m)

## Return a rotation matrix, given three angles in the order: XYZ
#
#  @see http://www.chrobotics.com/library/understanding-euler-angles
#  @param angles a list with angle x, y and z.
#  @return a matrix X * Y * Z
#
def rotateXYZ(angles):

	return rotateZYX(angles).T

## Returns a rotation matrix about a given axis.
#
#  @param angle rotation angle in degrees.
#  @param axis: 0 - x, 1 - y, 2 - z
#  @return rotation matrix.
#
def getRotationMatrix(angle, axis):
	if axis == 0:
		return rotate(angle, 1.0, 0.0, 0.0)
	elif axis == 1:
		return rotate(angle, 0.0, 1.0, 0.0)
	else:
		return rotate(angle, 0.0, 0.0, 1.0)

## Main program for testing.
def main():
    t = translate (3,4,5)

    print ("t = translation matrix =\n%s (type = %s)" % (t,type(t)) )

    r = rotate(90, 1,1,1)
    print ("r = rotation matrix =\n%s (type = %s)" % (r,type(r)) )

    m = t*r
    print ("t*r =\n%s (type = %s)" % (m,type(m)) )

if __name__=="__main__":
    sys.exit(main())
