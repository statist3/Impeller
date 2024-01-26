# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 12:33:27 2012

@author: 44792
"""

import math
import numpy

Beta1=30.0
Beta2=30.0
h=20.0
bthickness=2
b2=8.0
hb=((h-b2)/2.0)*0.9
nBlades=7
reverse=-1
r1=18.0
r2=70.0    #r2+1.5mm
extraRin=-2 #round((b2/1.0),2)
extraRout=0 #-5.
fileName='./impellerInlet_importSalome.py'


def logLine(startAngle):
    angle=startAngle
    Cord=[]
    angleList=[]
    beta=Beta1
    
    for i in numpy.arange(r1+bthickness,r2,0.5):
        
        if i==r1:
            pass
        elif i==r2:
            angle=(math.degrees((math.log(float(i)/float(i-0.5))/math.tan(math.radians(beta*reverse)))+math.radians(angle)))
            angleList.append(angle)
            x=float(i)*math.cos(math.radians(angle))
            y=float(i)*math.sin(math.radians(angle))
            Cord.append([x,y])
        else:
            angle=(math.degrees((math.log(float(i)/float(i-0.5))/math.tan(math.radians(beta*reverse)))+math.radians(angle)))
            beta+=(Beta2-Beta1)/(len(numpy.arange(r1,r2+0.5,0.5))-1)
            angleList.append(angle)
            x=float(i)*math.cos(math.radians(angle))
            y=float(i)*math.sin(math.radians(angle))
            Cord.append([x,y])
    return Cord

def writeCurve(l1,fileName,number):
    f = open(fileName, 'a')
    v=[]
    k=number
    
    for i,number in enumerate(l1):
        newVertex = ('V'+str(k)+'_'+str(i)+' = geompy.MakeVertex('+",".join(format(x, "10.5f") for x in number)+',0.0)\n')
        f.write(newVertex)
        if i>0:
            pass
        else:
            f.write('geompy.addToStudy(V'+str(k)+'_'+str(i)+', "V'+str(k)+'_'+str(i)+'" ) \n')
        v.append('V'+str(k)+'_'+str(i))
        
    f.write('Curve_'+str(k)+' = geompy.MakeInterpol(['+",".join(v)+'])\n')
    f.write('geompy.addToStudy(Curve_'+str(k)+', "Curve_'+str(k)+'" ) \n')
    f.close()
    
def writeStdSalome(fileName):
    f = open(fileName, 'w')
    f.write("""# -*- coding: utf-8 -*-
	
import sys
import salome
import os

fpath = os.path.dirname(sys.argv[0])

# INIT THE SALOME PART
salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, fpath)

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

geompy = geomBuilder.New()

# STL finenes
stlsize = 0.1


""")
    f.write('O = geompy.MakeVertex(0, 0, 0)\n')
    f.write('OX = geompy.MakeVectorDXDYDZ(1, 0, 0)\n')
    f.write('OY = geompy.MakeVectorDXDYDZ(0, 1, 0)\n')
    f.write('OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)\n')
    f.write('OZNeg = geompy.MakeVectorDXDYDZ(0, 0, -1)\n')
    f.write('geompy.addToStudy(O, "O" ) \n')
    f.write('geompy.addToStudy(OX, "OX" ) \n')
    f.write('geompy.addToStudy(OY, "OY" ) \n')
    f.write('geompy.addToStudy(OZ, "OZ" ) \n')
    f.write('geompy.addToStudy(OZNeg, "OZNeg" ) \n')
    f.close()
    
def writeArc(fileName,V1,V2,V3,number):
    f = open(fileName, 'a')
    f.write('Arc_'+str(number)+' = geompy.MakeArc('+V1+', '+V2+', '+V3+')\n')
    f.write('geompy.addToStudy(Arc_'+str(number)+', "Arc_'+str(number)+'" ) \n')
    f.close()

def writeFaceFromEdges(fileName,edgeList,fname):
    f = open(fileName, 'a')
    f.write(fname+' = geompy.MakeFaceWires(['+",".join(edgeList)+'], 1)\n')
    f.write('geompy.addToStudy('+fname+', "'+fname+'" ) \n')
    f.close()
    
def writeBase(filename,r2,b2):
    f = open(fileName, 'a')
    #f.write('shroudV = geompy.MakeVertex(0, 0, %d)\n' %(-h/2))
    f.write('hubV = geompy.MakeVertex(0, 0, %d)\n' %(h/2))
    f.write('inletV = geompy.MakeVertex(0, 0, %d)\n' %((h/2)+bthickness))
    #f.write('geompy.addToStudy(shroudV, "shroudV") \n')
    f.write('geompy.addToStudy(hubV, "hubV") \n')
    f.write('geompy.addToStudy(inletV, "inletV") \n')
    
    f.write('Cylinder_1 = geompy.MakeCylinder(hubV, OZ, %d, %d)\n' %(r1,bthickness))
    f.write('Cylinder_2 = geompy.MakeCylinder(inletV, OZ, %d, %d)\n' %(r1,r1*5))
    #f.write('Cylinder_3 = geompy.MakeCylinder(hubV, OZ, %d, %d)\n' %(r1,bthickness))
    #f.write('Cut_1 = geompy.MakeCutList(Cylinder_1, [Cylinder_3], True)\n')
    
    f.write('geompy.addToStudy( Cylinder_1, "Cylinder_1" )\n')
    f.write('geompy.addToStudy( Cylinder_2, "Cylinder_2" )\n')
    #f.write('geompy.addToStudy( Cylinder_3, "Cylinder_3" )\n')
    #f.write('geompy.addToStudy( Cut_1, "Cut_1" )\n')
    
    f.write('Extrusion_1 = geompy.MakePrismVecH2Ways(Curve_1, OZ, %d)\n' %(h/2))
    f.write('[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Extrusion_1, geompy.ShapeType["EDGE"], True)\n')
    f.write('[Vertex_1,Vertex_2,Vertex_3,Vertex_4] = geompy.ExtractShapes(Extrusion_1, geompy.ShapeType["VERTEX"], True)\n')
    f.write('Offset_1 = geompy.MakeOffset(Extrusion_1, %d)\n'% -bthickness)
    f.write('[Edge_5,Edge_6,Edge_7,Edge_8] = geompy.ExtractShapes(Offset_1, geompy.ShapeType["EDGE"], True)\n')
    f.write('[Vertex_5,Vertex_6,Vertex_7,Vertex_8] = geompy.ExtractShapes(Offset_1, geompy.ShapeType["VERTEX"], True)\n')
    f.write('Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_5)\n')
    f.write('Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_7)\n')

    f.write('Wire_2 = geompy.MakeWire([Line_1, Line_3, Edge_2, Edge_6], 1e-007)\n')
    f.write('geompy.addToStudy( Wire_2, "Wire_2" )\n')
    f.write('Fillet_1D_2 = geompy.MakeFillet1D(Wire_2, %s, [3, 4, 6, 8])\n' % ((bthickness/2.0)*0.99))
    f.write('Face_2 = geompy.MakeFaceWires([Fillet_1D_2], 1)\n')

    f.write('Multi_Rotation_1 = geompy.MultiRotate1DNbTimes(Face_2, None, %d)\n' % nBlades)
    
    f.write('geompy.addToStudyInFather( Extrusion_1, Edge_2, "Edge_2" )\n')
    f.write('geompy.addToStudyInFather( Extrusion_1, Edge_1, "Edge_1" )\n')
    f.write('geompy.addToStudyInFather( Extrusion_1, Edge_3, "Edge_3" )\n')
    f.write('geompy.addToStudyInFather( Extrusion_1, Edge_4, "Edge_4" )\n')
    f.write('geompy.addToStudyInFather( Offset_1, Edge_5, "Edge_5" )\n')
    f.write('geompy.addToStudyInFather( Offset_1, Edge_6, "Edge_6" )\n')
    f.write('geompy.addToStudyInFather( Offset_1, Edge_7, "Edge_7" )\n')
    f.write('geompy.addToStudyInFather( Offset_1, Edge_8, "Edge_8" )\n')
    f.write('geompy.addToStudyInFather( Extrusion_1, Vertex_1, "Vertex_1" )\n')
    f.write('geompy.addToStudyInFather( Extrusion_1, Vertex_2, "Vertex_2" )\n')
    f.write('geompy.addToStudyInFather( Extrusion_1, Vertex_3, "Vertex_3" )\n')
    f.write('geompy.addToStudyInFather( Extrusion_1, Vertex_4, "Vertex_4" )\n')
    f.write('geompy.addToStudyInFather( Offset_1, Vertex_5, "Vertex_5" )\n')
    f.write('geompy.addToStudyInFather( Offset_1, Vertex_6, "Vertex_6" )\n')
    f.write('geompy.addToStudyInFather( Offset_1, Vertex_7, "Vertex_7" )\n')
    f.write('geompy.addToStudyInFather( Offset_1, Vertex_8, "Vertex_8" )\n')
    f.write('geompy.addToStudy( Line_1, "Line_1" )\n')
    f.write('geompy.addToStudy( Line_3, "Line_3" )\n')
    f.write('geompy.addToStudy( Face_2, "Face_2" )\n')
    f.write('geompy.addToStudy( Multi_Rotation_1, "Multi-Rotation_1" )\n')
    
    
    f.write('LocalCS_1 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 0, 1)\n')
    f.write('sk = geompy.Sketcher2D()\n')
    f.write('sk.addPoint(0.000000, %d)\n' %(-h/2))
    f.write('sk.addSegmentAbsolute(%d, %d)\n' %(r1, -h/2))
    f.write('sk.addSegmentAbsolute(%d, %d)\n' %(r2, -h/2))
    f.write('sk.addSegmentAbsolute(%d, %d)\n' %(r2+5, -h/2))
    f.write('Sketch_1 = sk.wire(LocalCS_1)\n')
    f.write('Revolution_1 = geompy.MakeRevolution(Sketch_1, OZ, 360*math.pi/180.0)\n')
    
    f.write('geompy.addToStudy( LocalCS_1, "LocalCS_1" )\n')
    f.write('geompy.addToStudy( Sketch_1, "Sketch_1" )\n')
    f.write('geompy.addToStudy( Revolution_1, "Revolution_1" )\n')
    
    f.write('Cut_1 = geompy.MakeCutList(Revolution_1, [Multi_Rotation_1], True)\n')
    f.write('Extrusion_1 = geompy.MakePrismVecH(Cut_1, OZ, 20)\n')
    f.write('Fuse_1 = geompy.MakeFuseList([Cylinder_1, Extrusion_1], True, False)\n')
    
    f.write('geompy.addToStudy( Cut_1, "Cut_1" )\n')
    f.write('geompy.addToStudy( Extrusion_1, "Extrusion_1" )\n')
    f.write('geompy.addToStudy( Fuse_1, "Fuse_1" )\n')
    
    f.write('amiImpToInlet = geompy.CreateGroup(Fuse_1, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(amiImpToInlet, [10])\n')
    f.write('rotatingWallShroud = geompy.CreateGroup(Fuse_1, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(rotatingWallShroud, [12])\n')
    f.write('rotatingWallHub = geompy.CreateGroup(Fuse_1, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(rotatingWallHub, [377, 394])\n')
    f.write('wallImpeller = geompy.CreateGroup(Fuse_1, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(wallImpeller, [122, 3, 391])\n')
    f.write('rotatingWallImpeller = geompy.CreateGroup(Fuse_1, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(rotatingWallImpeller, [349, 299, 179, 144, 149, 279, 334, 244, 302, 274, 239, 169, 229, 162, 314, 319, 359, 209, 224, 194, 309, 289, 259, 189, 369, 329, 324, 249, 344, 204, 264, 337, 139, 154, 364, 219, 294, 254, 159, 174, 134, 127, 284, 197, 184, 214, 267, 232, 354])\n')
    f.write('amiImpToVolute = geompy.CreateGroup(Fuse_1, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(amiImpToVolute, [372])\n')
    
    f.write('geompy.addToStudyInFather( Fuse_1, amiImpToInlet, "amiImpToInlet" )\n')
    f.write('geompy.addToStudyInFather( Fuse_1, rotatingWallShroud, "rotatingWallShroud" )\n')
    f.write('geompy.addToStudyInFather( Fuse_1, rotatingWallHub, "rotatingWallHub" )\n')
    f.write('geompy.addToStudyInFather( Fuse_1, wallImpeller, "wallImpeller" )\n')
    f.write('geompy.addToStudyInFather( Fuse_1, rotatingWallImpeller, "rotatingWallImpeller" )\n')
    f.write('geompy.addToStudyInFather( Fuse_1, amiImpToVolute, "amiImpToVolute" )\n')
    
    f.write('amiInletToImp = geompy.CreateGroup(Cylinder_2, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(amiInletToImp, [12])\n')
    f.write('wallInlet = geompy.CreateGroup(Cylinder_2, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(wallInlet, [3])\n')
    f.write('inlet = geompy.CreateGroup(Cylinder_2, geompy.ShapeType["FACE"])\n')
    f.write('geompy.UnionIDs(inlet, [10])\n')
    
    f.write('geompy.addToStudyInFather( Cylinder_2, amiInletToImp, "amiInletToImp" )\n')
    f.write('geompy.addToStudyInFather( Cylinder_2, wallInlet, "wallInlet" )\n')
    f.write('geompy.addToStudyInFather( Cylinder_2, inlet, "inlet" )\n')

    
    f.close()

def writeFinal(filename):
    f = open(fileName, 'a')
    f.write('geompy.ExportSTL(amiImpToInlet, "Impeller/amiImpToInletMesh.stl", True, stlsize, False )\n')
    f.write('geompy.ExportSTL(rotatingWallShroud, "Impeller/rotatingWallShroudMesh.stl", True, stlsize, False )\n')
    f.write('geompy.ExportSTL(rotatingWallHub, "Impeller/rotatingWallHubMesh.stl", True, stlsize, False )\n')
    f.write('geompy.ExportSTL(wallImpeller, "Impeller/wallImpellerMesh.stl", True, stlsize, False )\n')
    f.write('geompy.ExportSTL(rotatingWallImpeller, "Impeller/rotatingWallImpellerMesh.stl", True, stlsize, False )\n')
    f.write('geompy.ExportSTL(amiImpToVolute, "Impeller/amiImpToVoluteMesh.stl", True, stlsize, False )\n')
    
    f.write('geompy.ExportSTL(amiInletToImp, "Inlet/amiInletToImpMesh.stl", True, stlsize, False )\n')
    f.write('geompy.ExportSTL(wallInlet, "Inlet/wallInletMesh.stl", True, stlsize, False )\n')
    f.write('geompy.ExportSTL(inlet, "Inlet/inletMesh.stl", True, stlsize, False )\n')
    f.close()

l1=numpy.array(logLine(0))

writeStdSalome(fileName)

writeCurve(l1,fileName,1)

writeBase(fileName,r2,b2)

writeFinal(fileName)


