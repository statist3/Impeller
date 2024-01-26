# -*- coding: utf-8 -*-

import numpy
import math
### Std Input###

filename='./Volute_importSalome.py'
D2=140+10           #Impeller diameter
DIN=50                #Outlet diameter
tongueThickness=8.0     #Tongue thickness
zMin=-10             #Volute Height -z
zMax=10              #Volute Height +x
outFactor=0.95          #% of outletarea that the tonguearea should be should be between 0.8-1
         

###Advanced Input###


tonguePush=0.5          #Push the tongue outwards this amount of mm
tongueGrowthConst=1.025  #Growth constant for the first 30 degree of volute, to minimise protrusion into the outlet pipe
extraTongue=2          #Extra tongue dimensions to make the flat section of the volute longer

outletX=1              #Extra X length to make the pipe not intersect with the volute
outletY=55             #Extra Y length to make the pipe not intersect with the volute
outletZ=0               #Pull the outlet in the z-direction
zNeut=0


###End of user input###
outletArea=int(math.pi*math.pow(DIN/2,2))
#print "outletArea %d mm^2"%outletArea

GrowthConst=1.02  # start guess of the growthConstant
Radius=D2/2
overHungGrow=1.0

angleList=[0,10,20,30,60,90,120,150,180,210,240,270,300,330,360]
lenAngleList=len(angleList)-1


def findCircleFit(p1x, p1y, p2x, p2y, p3x, p3y):
    """Function to find the circle parameters from three points"""
    alpha=numpy.linalg.det([[p1x,p1y,1],[p2x,p2y,1],[p3x,p3y,1]])
#    print alpha
    beta=numpy.linalg.det([[p1x**2+p1y**2,p1y,1],[p2x**2+p2y**2,p2y,1],[p3x**2+p3y**2,p3y,1]])
#    print beta
    gamma=numpy.linalg.det([[p1x**2+p1y**2,p1x,1],[p2x**2+p2y**2,p2x,1],[p3x**2+p3y**2,p3x,1]])
#    print gamma
    sigma=numpy.linalg.det([[p1x**2+p1y**2,p1x,p1y],[p2x**2+p2y**2,p2x,p2y],[p3x**2+p3y**2,p3x,p3y]])
#    print sigma
    
    if alpha==0:
        print("points are collinear, alpha set to 10^-30")
        alpha=10^-30
    
    x0=beta/(2.0*alpha)
#    print x0
    y0=-1.0*gamma/(2.0*alpha)
#    print y0
    r0=math.sqrt(x0**2+y0**2+(sigma/alpha))
    return x0,y0,r0
    
p1=[tongueThickness+Radius*math.pow(overHungGrow,lenAngleList+1)*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zMax]
p2=[tongueThickness+Radius*math.pow(GrowthConst,lenAngleList+1)*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zNeut]
p3=[tongueThickness+Radius*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zMin]

#findCircleFit(66.0,1,89.8,0,53.0,0)
#print p1[0],p1[1],p2[0],p2[1],p3[0],p3[1]
xT,zT,rT=findCircleFit(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
#print xT,zT,rT

tongueArea=math.pi*math.pow(rT,2)
#print "tongueArea %d mm^2"%tongueArea

#A simple goal seek that finds the optimal value of the GrowthConst
while tongueArea <= (outletArea*outFactor*0.98) or tongueArea >= (outletArea*outFactor*1.02):
    if GrowthConst<1.006642:
        GrowthConst=1.0006642
        break
    elif tongueArea <= outletArea*outFactor:
        #print "Smaller"
        GrowthConst=GrowthConst+0.000001
        #print "GrowthConst= "+str(GrowthConst)
        p1=[tongueThickness+Radius*math.pow(overHungGrow,lenAngleList+1)*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zMax]
        p2=[tongueThickness+Radius*math.pow(GrowthConst,lenAngleList+1)*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zNeut]
        p3=[tongueThickness+Radius*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zMin]
        xT,zT,rT=findCircleFit(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
        tongueArea=math.pi*math.pow(rT,2)
        #print "tongueArea {0} mm^2, targetArea {1} mm^2, outletArea {2} mm^2 \n".format(tongueArea,float(outletArea*outFactor),outletArea)
    elif tongueArea > outletArea*outFactor:
        #print "Larger"
        GrowthConst=GrowthConst-0.000001
        #print "GrowthConst= "+str(GrowthConst)
        p1=[tongueThickness+Radius*math.pow(overHungGrow,lenAngleList+1)*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zMax]
        p2=[tongueThickness+Radius*math.pow(GrowthConst,lenAngleList+1)*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zNeut]
        p3=[tongueThickness+Radius*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zMin]
        xT,zT,rT=findCircleFit(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
        tongueArea=math.pi*math.pow(rT,2)
        #print "tongueArea {0} mm^2, targetArea {1} mm^2, outletArea {2} mm^2 \n".format(tongueArea,float(outletArea*outFactor),outletArea)
    else:
        #print "Match"
        GrowthConst=GrowthConst
        #print 
        p1=[tongueThickness+Radius*math.pow(overHungGrow,lenAngleList+1)*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zMax]
        p2=[tongueThickness+Radius*math.pow(GrowthConst,lenAngleList+1)*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zNeut]
        p3=[tongueThickness+Radius*math.cos(math.radians(angleList[lenAngleList]))+extraTongue*math.cos(math.radians(angleList[lenAngleList])),zMin]
        xT,zT,rT=findCircleFit(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
        tongueArea=math.pi*math.pow(rT,2)
        #print "tongueArea {0} mm^2, targetArea {1} mm^2, outletArea {2} mm^2 \n".format(tongueArea,float(outletArea*outFactor),outletArea)

        

print("GrowthConst= "+str(GrowthConst))

midSpan=[]
topSpan=[]
topSpanCirc=[]
botSpan=[]
tongueVertex=[]
arcList=[]
tongueArcList=[]

for i,j in enumerate(angleList):
    i+=1
    if j==360:
        midSpan.append([tongueThickness+Radius*math.pow(GrowthConst,i)*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.pow(GrowthConst,i)*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zNeut])
        topSpan.append([tongueThickness+Radius*math.pow(overHungGrow,i)*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMax])
        topSpanCirc.append([tongueThickness+Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMax])
        botSpan.append([tongueThickness+Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMin])
        tongueVertex.append([Radius*math.cos(math.radians(j))+tongueThickness*math.cos(math.radians(j))+(extraTongue-1)*math.cos(math.radians(j)),Radius*math.pow(GrowthConst,1)*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zNeut]) 

    elif j==0:
        midSpan.append([Radius+tonguePush*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.pow(GrowthConst,i)*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zNeut])
        topSpan.append([Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMax])
        topSpanCirc.append([Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMax])
        botSpan.append([Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMin])
    elif j<=30:
        midSpan.append([Radius*math.pow(tongueGrowthConst,i)*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.pow(tongueGrowthConst,i)*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zNeut])
        topSpan.append([Radius*math.pow(overHungGrow,i)*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.pow(overHungGrow,i)*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMax])
        topSpanCirc.append([Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMax])
        botSpan.append([Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMin])
    else:
        midSpan.append([Radius*math.pow(GrowthConst,i)*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.pow(GrowthConst,i)*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zNeut])
        topSpan.append([Radius*math.pow(overHungGrow,i)*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.pow(overHungGrow,i)*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMax])
        topSpanCirc.append([Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMax])
        botSpan.append([Radius*math.cos(math.radians(j))+extraTongue*math.cos(math.radians(j)),Radius*math.sin(math.radians(j))+extraTongue*math.sin(math.radians(j)),zMin])

f = open(filename, 'w')
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
f.write('geompy.addToStudy(O, "O" ) \n')
f.write('geompy.addToStudy(OX, "OX" ) \n')
f.write('geompy.addToStudy(OY, "OY" ) \n')
f.write('geompy.addToStudy(OZ, "OZ" ) \n')

for i,number in enumerate(tongueVertex):
    newVertex = ('V'+str(i)+'Tongue = geompy.MakeVertex('+",".join(format(x, "10.5f") for x in number)+')\n')
#    print newVertex
    f.write(newVertex)
    #f.write('geompy.addToStudy(V'+str(i)+'Tongue, "V'+str(i)+'Tongue" ) \n')

for i,number in enumerate(midSpan):
    newVertex = ('V'+str(i)+'Mid = geompy.MakeVertex('+",".join(format(x, "10.5f") for x in number)+')\n')
#    print newVertex
    f.write(newVertex)
#    f.write('geompy.addToStudy(V'+str(i)+'Mid, "V'+str(i)+'Mid" ) \n')
    
for i,number in enumerate(topSpan):
    newVertex = ('V'+str(i)+'Top = geompy.MakeVertex('+",".join(format(x, "10.5f") for x in number)+')\n')
#    print newVertex
    f.write(newVertex)
#    f.write('geompy.addToStudy(V'+str(i)+'Top, "V'+str(i)+'Top" ) \n')
    
for i,number in enumerate(topSpanCirc):
    newVertex = ('V'+str(i)+'TopCirc = geompy.MakeVertex('+",".join(format(x, "10.5f") for x in number)+')\n')
#    print newVertex
    f.write(newVertex)
#    f.write('geompy.addToStudy(V'+str(i)+'TopCirc, "V'+str(i)+'TopCirc" ) \n')
    
for i,number in enumerate(botSpan):
    newVertex = ('V'+str(i)+'Bot = geompy.MakeVertex('+",".join(format(x, "10.5f") for x in number)+')\n')
#    print newVertex
    f.write(newVertex)
#    f.write('geompy.addToStudy(V'+str(i)+'Bot, "V'+str(i)+'Bot" ) \n')
    
for i,j in enumerate(angleList):
    if j==360:
        newArc = ('Arc'+str(i)+' = geompy.MakeArc(V'+str(i)+'Top, V'+str(i)+'Mid, V'+str(i)+'Bot)\n')
        arcList.append("Arc"+str(i))
        f.write(newArc)
        #f.write('geompy.addToStudy(Arc'+str(i)+', "Arc'+str(i)+'" ) \n')
        
        newArc = ('ArcTongue = geompy.MakeArc(V'+str(i)+'TopCirc, V0Tongue, V'+str(i)+'Bot)\n')
        f.write(newArc)
        #f.write('geompy.addToStudy(ArcTongue, "ArcTongue" ) \n')
    elif j<30:
        newArc = ('Arc'+str(i)+' = geompy.MakeArc(V'+str(i)+'Top, V'+str(i)+'Mid, V'+str(i)+'Bot)\n')
        tongueArcList.append("Arc"+str(i))
        f.write(newArc)
        #f.write('geompy.addToStudy(Arc'+str(i)+', "Arc'+str(i)+'" ) \n')
    elif j==30:
        newArc = ('Arc'+str(i)+' = geompy.MakeArc(V'+str(i)+'Top, V'+str(i)+'Mid, V'+str(i)+'Bot)\n')
        tongueArcList.append("Arc"+str(i))
        arcList.append("Arc"+str(i))
        f.write(newArc)
        #f.write('geompy.addToStudy(Arc'+str(i)+', "Arc'+str(i)+'" ) \n')
    else:
        newArc = ('Arc'+str(i)+' = geompy.MakeArc(V'+str(i)+'Top, V'+str(i)+'Mid, V'+str(i)+'Bot)\n')
        arcList.append("Arc"+str(i))
        f.write(newArc)
        #f.write('geompy.addToStudy(Arc'+str(i)+', "Arc'+str(i)+'" ) \n')
    

f.write('VoluteTongueWall=geompy.MakeCompound(['+",".join(tongueArcList)+'])\n')
f.write('VoluteWall=geompy.MakeCompound(['+",".join(arcList)+'])\n')
#f.write('geompy.addToStudy(VoluteWall, "VoluteWall" ) \n')

f.write('TongueSweep = geompy.MakeFilling(VoluteTongueWall, 9, 10, 0.0001, 0.0001, 0, GEOM.FOM_AutoCorrect)\n')
f.write('geompy.addToStudy(TongueSweep, "TongueSweep" ) \n')

f.write('VoluteSweep = geompy.MakeFilling(VoluteWall, 9, 10, 0.0001, 0.0001, 0, GEOM.FOM_AutoCorrect)\n')
f.write('geompy.addToStudy(VoluteSweep, "VoluteSweep" ) \n')

f.write('Line_1 = geompy.MakeLineTwoPnt(V0Bot, V'+str(lenAngleList)+'Bot)\n')
f.write('Line_2 = geompy.MakeLineTwoPnt(V0TopCirc, V'+str(lenAngleList)+'TopCirc)\n')
#f.write('Line_3 = geompy.MakeLineTwoPnt(V'+str(lenAngleList)+'TopCirc, V'+str(lenAngleList)+'Top)\n')
#f.write('Line_4 = geompy.MakeLineTwoPnt(V0Bot, V0TopCirc)\n')
f.write('geompy.addToStudy(Line_1, "Line_1" ) \n')
f.write('geompy.addToStudy(Line_2, "Line_2" ) \n')
#f.write('geompy.addToStudy(Line_3, "Line_3" ) \n')
#f.write('geompy.addToStudy(Line_4, "Line_4" ) \n')

f.write('[Edge_1,Edge_2] = geompy.SubShapes(TongueSweep, [8, 3]) \n')

f.write('[Edge_3,Edge_4] = geompy.SubShapes(VoluteSweep, [8, 3]) \n')

f.write('Vertex_1 = geompy.MakeVertex(0, 0, '+str(zMax)+')\n')
f.write('Circle_1 = geompy.MakeCircle(Vertex_1, None, '+str(Radius)+')\n')
f.write('shroudFace = geompy.MakeFaceWires([Line_2, Circle_1, Edge_2,Edge_4], 0)\n')
f.write('Vertex_2 = geompy.MakeVertex(0, 0, '+str(zMin)+')\n')
f.write('Circle_2 = geompy.MakeCircle(Vertex_2, None, '+str(Radius)+')\n')
f.write('hubFace = geompy.MakeFaceWires([Line_1, Circle_2, Edge_1,Edge_3], 0)\n')

f.write('geompy.addToStudy(shroudFace, "shroudFace" ) \n')
f.write('geompy.addToStudy(hubFace, "hubFace" ) \n')

f.write('toungeCenterVertex = geompy.MakeVertex('+str(xT+tongueThickness)+', 0, '+str(zT)+')\n')
#f.write('geompy.addToStudy(toungeCenterVertex, "toungeCenterVertex" ) \n')

f.write('outletVertex = geompy.MakeVertexWithRef(toungeCenterVertex, '+str(outletX)+', '+str(outletY)+', '+str(outletZ)+')\n')
#f.write('geompy.addToStudy(outletVertex, "outletVertex" ) \n')

f.write('yVector = geompy.MakeLineTwoPnt(toungeCenterVertex, outletVertex)\n')
f.write('Circle_3 = geompy.MakeCircle(outletVertex, yVector, '+str(DIN/2)+')\n')
#f.write('geompy.addToStudy(Circle_3, "Circle_3" ) \n')

f.write('TongueFace = geompy.MakeFaceWires([Line_1, Line_2, Arc0, ArcTongue], 1)\n')
f.write('geompy.addToStudy(TongueFace, "TongueFace" ) \n')

x1=xT+tongueThickness
y1=0
x2=xT+tongueThickness+outletX
y2=outletY

m=(y2-y1)/(x2-x1)    

x=x2    
y=m*(x-x2)+y2
angleC = math.degrees(math.atan(x/y))
atan = math.degrees(math.atan(outletY/outletX))
angle=abs(((90-atan)-angleC))

if outletX>0:
    f.write('outVert1 = geompy.MakeVertexOnCurve(Circle_3, 0.6)\n')
    f.write('Extrusion_1 = geompy.MakePrismDXDYDZ(outVert1, 0, 0, '+str(DIN*-1.0)+')\n')
    f.write('Translation_1 = geompy.MakeTranslation(Extrusion_1, 0, 0, -1)\n')
    f.write('outVert2 = geompy.MakeVertexOnLinesIntersection(Circle_3, Translation_1)\n')
    f.write('outArc1=geompy.MakeArcCenter(outletVertex, outVert1, outVert2,False)\n')
    f.write('outArc2=geompy.MakeArcCenter(outletVertex, outVert1, outVert2,True)\n')
#        f.write('geompy.addToStudy(outVert1, "outVert1" ) \n')
#        f.write('geompy.addToStudy(outVert2, "outVert2" ) \n')
#        f.write('geompy.addToStudy(outArc1, "outArc1" ) \n')
#        f.write('geompy.addToStudy(outArc2, "outArc2" ) \n') 
elif outletX==0:
    f.write('outVert1 = geompy.MakeVertexOnCurve(Circle_3, 0.75)\n')
    f.write('Extrusion_1 = geompy.MakePrismDXDYDZ(outVert1, 0, 0, '+str(DIN*-1.0)+')\n')
    f.write('Translation_1 = geompy.MakeTranslation(Extrusion_1, 0, 0, -1)\n')
    f.write('outVert2 = geompy.MakeVertexOnLinesIntersection(Circle_3, Translation_1)\n')
    f.write('outArc1=geompy.MakeArcCenter(outletVertex, outVert1, outVert2,False)\n')
    f.write('outArc2=geompy.MakeArcCenter(outletVertex, outVert1, outVert2,True)\n')
#        f.write('geompy.addToStudy(outVert1, "outVert1" ) \n')
#        f.write('geompy.addToStudy(outVert2, "outVert2" ) \n')
#        f.write('geompy.addToStudy(outArc1, "outArc1" ) \n')
#        f.write('geompy.addToStudy(outArc2, "outArc2" ) \n') 

f.write('outerTongueComp = geompy.MakeCompound([Arc'+str(lenAngleList)+', outArc2])\n')
f.write('innerTongueComp = geompy.MakeCompound([ArcTongue, outArc1])\n')
f.write('outerTongueFace = geompy.MakeFilling(outerTongueComp, 9, 10, 0.0001, 0.0001, 0, GEOM.FOM_AutoCorrect)\n')
f.write('innerTongueFace = geompy.MakeFilling(innerTongueComp, 9, 10, 0.0001, 0.0001, 0, GEOM.FOM_AutoCorrect)\n')

f.write('geompy.addToStudy(outerTongueFace, "outerTongueFace" ) \n')
f.write('geompy.addToStudy(innerTongueFace, "innerTongueFace" ) \n')

f.write('Wire_1 = geompy.MakeWire([outArc1, outArc2], 1e-07)\n')
f.write('geompy.addToStudy(Wire_1, "Wire_1" ) \n')

f.write('Wire_2 = geompy.MakeSketcher("Sketcher:F '+str(x2)+' '+str(y2)+':T '+str(outletX/2.0)+' '+str(outletY/2.0)+':D '+str(x+outletX/2.0)+' '+str(y+outletY/2.0)+':L '+str(DIN*10.0)+'", [0, 0, 0, 0, 0, 1, 1, 0, -0])\n')
f.write('geompy.addToStudy(Wire_2, "Wire_2" ) \n')
#    f.write('Wire_2 = geompy.MakeSketcher("Sketcher:F '+str(x2)+' '+str(y2)+':D '+str(outletX)+' '+str(outletX)+':L '+str(DIN)+':R 0:C '+str(DIN*-1.0)+' '+str(angle)+':L '+str(DIN*3.0)+'", [0, 0, 0, 0, 0, 1, 1, 0, -0])\n')
f.write('Fillet_1D_1 = geompy.MakeFillet1D(Wire_2, '+str(DIN)+', [4])\n')
#    "Sketcher:F 102.953512 20.000000:D 1.000000 20.000000:L 20.000000:R 0:C -20.000000 66.091000:L 200.000000"

#    f.write('Wire_2 = geompy.MakeSketcher("Sketcher:F '+str(xT+tongueThickness)+' '+str(outletY)+':TT '+str(xT+tongueThickness)+' '+str(outletY+1.0)+':AA '+str(xT+tongueThickness+2.0)+' '+str(xT+tongueThickness+2.0)+':TT '+str(DIN*4.0)+' '+str(DIN*4.0)+'", [0, 0, 0, 0, 0, 1, 1, 0, -0])\n')
f.write('geompy.addToStudy(Wire_2, "Wire_2" ) \n')
f.write('geompy.addToStudy(Fillet_1D_1, "Fillet_1D_1" ) \n')
f.write('Pipe_1 = geompy.MakePipe(Wire_1, Fillet_1D_1)\n')
f.write('geompy.addToStudy(Pipe_1, "Pipe_1" ) \n')      
    
f.write('Shell = geompy.MakeShell([TongueSweep,VoluteSweep, shroudFace, hubFace, TongueFace, outerTongueFace, innerTongueFace,Pipe_1]) \n')
f.write('geompy.addToStudy(Shell, "Shell" ) \n')
f.write('Vertex_1 = geompy.MakeVertex('+str(x+outletX/2.0)+', '+str(y+outletY/2.0)+', 0)\n')
f.write('Vertex_2 = geompy.MakeVertex(0, 0, 0)\n')
f.write('Vertex_3 = geompy.MakeVertex('+str(x+outletX/2.0+(DIN*10))+', 0, 0)\n')

f.write('geompy.addToStudy(Vertex_1, "Vertex_1" ) \n')
f.write('geompy.addToStudy(Vertex_2, "Vertex_2" ) \n')
f.write('refineTongue = geompy.TrsfOp.RotateThreePoints(TongueSweep, Vertex_2, Vertex_1, toungeCenterVertex)\n')
f.write('geompy.TrsfOp.RotateThreePoints(Shell, Vertex_2, Vertex_1, toungeCenterVertex)\n')



f.write('Vertex_4 = geompy.MakeVertex(%d, 0, %d)\n'%(D2/2,zMin))
f.write('Vertex_5 = geompy.MakeVertex(%d, 0, %d)\n'%(D2/2,zMax))
f.write('Line_5 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_5)\n')
f.write('Revolution_1 = geompy.MakeRevolution(Line_5, OZ, 360*math.pi/180.0)\n')
f.write('geompy.addToStudy(Revolution_1, "Revolution_1" ) \n')

f.write('Shell_2 = geompy.MakeShell([Shell,Revolution_1]) \n')
f.write('geompy.addToStudy(Shell_2, "Shell_2" ) \n')

f.write('Vertex_6 = geompy.MakeVertex(%d, 0, 0)\n'%(DIN*10))
f.write('Cylinder_1 = geompy.MakeCylinder(Vertex_6, OX, %d, 300)\n' %(DIN*2))
f.write('Cut_1 = geompy.MakeCutList(Shell_2, [Cylinder_1], True)\n')
f.write('SuppressHoles_1 = geompy.SuppressHoles(Cut_1, [])\n')

f.write('geompy.addToStudy(SuppressHoles_1, "SuppressHoles_1" ) \n')

f.write('amiVoluteToImp = geompy.CreateGroup(SuppressHoles_1, geompy.ShapeType["FACE"])\n')
f.write('geompy.UnionIDs(amiVoluteToImp, [69])\n')
f.write('wallTongue = geompy.CreateGroup(SuppressHoles_1, geompy.ShapeType["FACE"])\n')
f.write('geompy.UnionIDs(wallTongue, [43, 36, 29])\n')
f.write('wallTopBottom = geompy.CreateGroup(SuppressHoles_1, geompy.ShapeType["FACE"])\n')
f.write('geompy.UnionIDs(wallTopBottom, [50, 58])\n')
f.write('wallVolute = geompy.CreateGroup(SuppressHoles_1, geompy.ShapeType["FACE"])\n')
f.write('geompy.UnionIDs(wallVolute, [72, 66])\n')
f.write('wallOutlet = geompy.CreateGroup(SuppressHoles_1, geompy.ShapeType["FACE"])\n')
f.write('geompy.UnionIDs(wallOutlet, [78, 75, 81, 22, 8, 15])\n')
f.write('outlet = geompy.CreateGroup(SuppressHoles_1, geompy.ShapeType["FACE"])\n')
f.write('geompy.UnionIDs(outlet, [2])\n')

f.write('geompy.addToStudyInFather( SuppressHoles_1, amiVoluteToImp, "amiVoluteToImp" )\n')
f.write('geompy.addToStudyInFather( SuppressHoles_1, wallTongue, "wallTongue" )\n')
f.write('geompy.addToStudyInFather( SuppressHoles_1, wallTopBottom, "wallTopBottom" )\n')
f.write('geompy.addToStudyInFather( SuppressHoles_1, wallVolute, "wallVolute" )\n')
f.write('geompy.addToStudyInFather( SuppressHoles_1, wallOutlet, "wallOutlet" )\n')
f.write('geompy.addToStudyInFather( SuppressHoles_1, outlet, "outlet" )\n')


f.write('geompy.ExportSTL(amiVoluteToImp, "Volute/amiVoluteToImpMesh.stl", True, stlsize, False )\n')
f.write('geompy.ExportSTL(wallTongue, "Volute/wallTongueMesh.stl", True, stlsize, False )\n')
f.write('geompy.ExportSTL(wallTopBottom, "Volute/wallTopBottomMesh.stl", True, stlsize, False )\n')
f.write('geompy.ExportSTL(wallVolute, "Volute/wallVoluteMesh.stl", True, stlsize, False )\n')
f.write('geompy.ExportSTL(wallOutlet, "Volute/wallOutletMesh.stl", True, stlsize, False )\n')
f.write('geompy.ExportSTL(outlet, "Volute/outletMesh.stl", True, stlsize, False )\n')


f.close()
