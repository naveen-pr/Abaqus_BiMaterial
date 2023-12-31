"""
Written by Naveen Prakash
"""

# -*- coding: mbcs -*-
import numpy as np
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqus import *
from abaqusConstants import *
from odbAccess import *

## Inputs Params - MPa and mm

length=50
mat1_thk=0.5
mat1_E=70e3
mat1_v=0.22
mat1_cte=5e-6
#
mat2_thk=0.9
mat2_E=105000.0
mat2_v=0.22
mat2_cte=2.877175000000256e-07
#
delta_temp=100
app_stress=0.01
#
num_elem_length=20
num_elem_thk=2 #per layer
#
eps=length*1e-3
#
first_order=False
enhanced=False
#
jobname='BiMaterial'

## Base geometry -----

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# My model alias
myModel = mdb.models['Model-1'] 

mySketch = myModel.ConstrainedSketch(name='__profile__', sheetSize=100.0)
mySketch.rectangle(point1=(0.0, 0.0), point2=(length, (mat1_thk + mat2_thk) ))

mySketch.vertices.findAt((0.0, 0.0))
mySketch.FixedConstraint(entity=mySketch.vertices.findAt((0.0, 0.0), ))

myModel.Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=DEFORMABLE_BODY)
myModel.parts['Part-1'].BaseShell(sketch=mySketch)
del mySketch

# My part alias
myPart = myModel.parts['Part-1']

## Partition into two materials

mySketch = myModel.ConstrainedSketch(gridSpacing=1.0, name='__profile__', 
			sheetSize=40.11, transform=
			myPart.MakeSketchTransform(
			sketchPlane=myPart.faces.findAt((eps, eps, 0.0), ), sketchPlaneSide=SIDE1, 
			sketchOrientation=RIGHT, origin=(0., 0., 0.)))		
myPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mySketch)
mySketch.Line(point1=(0.0, mat2_thk), point2=(length, mat2_thk))
myPart.PartitionFaceBySketch(faces=myPart.faces.findAt(((eps, eps, 0.0), )), sketch=mySketch)
del mySketch

## Materials ------

myModel.Material(name='Mat1')
myModel.materials['Mat1'].Elastic(table=((mat1_E, mat1_v), ))
myModel.materials['Mat1'].Expansion(table=((mat1_cte, ), ))

myModel.Material(name='Mat2')
myModel.materials['Mat2'].Elastic(table=((mat2_E, mat2_v), ))
myModel.materials['Mat2'].Expansion(table=((mat2_cte, ), ))

## Create Sets -----

# Create regions
def createRegion(*args):
    
    if args[0] == 'face' or args[0] == 'faces':
        region = myPart.faces.findAt((args[1], ),)
        for coord in args[1:]:
            region = region + myPart.faces.findAt((coord, ), )
            print(coord)
    
    elif args[0] == 'edge' or args[0] == 'edges':
        region = myPart.edges.findAt((args[1], ),)
        for coord in args[1:]:
            region = region + myPart.edges.findAt((coord, ), )
    
    elif args[0] == 'vertex' or args[0] == 'vertices':
        region = myPart.vertices.findAt((args[1], ),)
        for coord in args[1:]:
            region = region + myPart.vertices.findAt((coord, ), )
        
    return region

reg_Mat2        = createRegion('faces', (length/2., mat2_thk/2., 0.0) )
reg_Mat1        = createRegion('faces', (length/2., mat2_thk + mat1_thk/2., 0.0) )
reg_All         = createRegion('faces', (length/2., mat2_thk/2., 0.0), (length/2., mat2_thk + mat1_thk/2., 0.0) )
reg_XSymm       = createRegion('edges', (0.0, mat2_thk/2., 0.0), (0.0, mat2_thk + mat1_thk/2., 0.0) )
reg_TopSurf     = createRegion('edges', (length/2., mat1_thk + mat2_thk, 0.0) )
reg_BotSurf     = createRegion('edges', (length/2., 0.0, 0.0) )
reg_Interface   = createRegion('edges', (length/2., mat2_thk, 0.0) )
reg_Tip         = createRegion('vertices', (length, 0.0, 0.0) )
reg_SuppCorner  = createRegion('vertices', (0.0, mat2_thk + mat1_thk, 0.0) )
reg_Encastre    = createRegion('vertices', (0.0, 0.0, 0.0) )

# Create sets
myPart.Set(faces=reg_Mat2, name='Mat2')
myPart.Set(faces=reg_Mat1, name='Mat1')
myPart.Set(edges=reg_XSymm, name='XSymm')
myPart.Set(vertices=reg_Encastre, name='Encastre')
myPart.Set(vertices=reg_Tip, name='Tip')
myPart.Set(vertices=reg_SuppCorner, name='SupportCorner')
myPart.Set(edges=reg_TopSurf, name='Top Surface')
myPart.Set(edges=reg_BotSurf, name='Bottom Surface')
myPart.Set(edges=reg_Interface, name='Interface')
myPart.Set(faces=reg_Mat1+reg_Mat2, name='All')
myPart.Surface(side1Edges=reg_TopSurf, name='Top Surface')

## Assign Sections -----
myModel.HomogeneousSolidSection(material='Mat1', name='Section-1', thickness=1.0)
myModel.HomogeneousSolidSection(material='Mat2', name='Section-2', thickness=1.0)

myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
    myPart.sets['Mat1'], sectionName='Section-1'
    , thicknessAssignment=FROM_SECTION)
myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
    myPart.sets['Mat2'], sectionName='Section-2'
    , thicknessAssignment=FROM_SECTION)
	
## Create Instance -----
myModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
myInstance = myModel.rootAssembly.Instance(dependent=ON, name='Part-1-1', part=myPart)

## Steps -----

# Initial step
myModel.Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(0.0, ), name='Temperature', region=myInstance.sets['All'])
	
myModel.XsymmBC(createStepName='Initial', localCsys=None, name=
    'XSymm', region=myInstance.sets['XSymm'])
mdb.models['Model-1'].EncastreBC(createStepName='Initial', localCsys=None, 
    name='Encastre', region=myInstance.sets['Encastre'])

# Static step
myModel.StaticStep(name='Step-1', nlgeom=ON, previous='Initial')
myModel.predefinedFields['Temperature'].setValuesInStep(magnitudes=(delta_temp, ), stepName='Step-1')
myModel.Pressure(amplitude=UNSET, createStepName='Step-1', 
	distributionType=UNIFORM, field='', magnitude=app_stress, name='Load-1', region=myInstance.surfaces['Top Surface'])

## Mesh -----

myPart.seedEdgeByNumber(constraint=FINER, edges=reg_TopSurf, number=num_elem_length)
myPart.seedEdgeByNumber(constraint=FINER, edges=reg_BotSurf, number=num_elem_length)
myPart.seedEdgeByNumber(constraint=FINER, edges=reg_Interface, number=num_elem_length)
myPart.seedEdgeByNumber(constraint=FINER, edges=reg_XSymm, number=num_elem_thk)

if first_order:
    myPart.setElementType(elemTypes=(ElemType(
	elemCode=CPS4R, elemLibrary=STANDARD), ElemType(
	elemCode=CPS3, elemLibrary=STANDARD)), regions=(reg_All, )) 
    
    if enhanced:
        myPart.setElementType(elemTypes=(ElemType(
   		elemCode=CPS4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
   		hourglassControl=ENHANCED, distortionControl=DEFAULT), ElemType(
   		elemCode=CPS3, elemLibrary=STANDARD)), regions=(reg_All, ))
else:
    myPart.setElementType(elemTypes=(ElemType(
    elemCode=CPS8R, elemLibrary=STANDARD), ElemType(elemCode=CPS6M, 
    elemLibrary=STANDARD)), regions=(reg_All, )) 

myPart.generateMesh()

## History Outputs -----
myModel.HistoryOutputRequest(createStepName='Step-1', frequency=
    LAST_INCREMENT, name='TipDisplacement', rebar=EXCLUDE, region=myInstance.sets['Tip'], 
    sectionPoints=DEFAULT, variables=('U2', ))

# Create job
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, name=
    jobname, nodalOutputPrecision=SINGLE, queue=None, resultsFormat=ODB, 
    scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
myModel.rootAssembly.regenerate()
mdb.saveAs(jobname)

job = mdb.jobs[jobname]
job.writeInput()
job.submit(consistencyChecking=OFF)
job.waitForCompletion()

#Open odb file
session.openOdb(name=jobname + '.odb', readOnly=True)
odb = session.odbs[jobname + '.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)

data1 = session.xyDataListFromField(odb=odb, outputPosition=NODAL, 
			variable=(('U', NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=("PART-1-1.TIP", ))

data2 = session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('S', 
        INTEGRATION_POINT, ((COMPONENT, 'S11'), )), ), nodeSets=(
        "PART-1-1.SUPPORTCORNER", ))

tip_disp    = data1[0][-1][-1]			
max_stress  = data2[0][-1][-1]

file_path = 'results.dat'

# The line you want to append
new_line = str(mat2_thk) + '\t' + str(mat2_E) + '\t' + str(mat2_cte) + '\t' + str(tip_disp) + '\t' + str(max_stress)

# Open the file in append mode and write the new line
with open(file_path, 'a') as file:
    file.write(new_line + '\n')

