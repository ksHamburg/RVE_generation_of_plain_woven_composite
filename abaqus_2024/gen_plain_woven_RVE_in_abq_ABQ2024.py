''' This abaqus-python program creates a fully periodic (geometry and mesh) RVE 
of plain woven composite.
Materialorientation needs to be implemented into the model via GUI afterwards.
If error messages occur during the change of geometry parameter modification 
overlap might occur.

If error messages pop up after geometry creation completed, they are most likely 
due to meshing. May be adjusting the seed density helps.
'''

try:
    session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
except NameError:
    print("No session available")
    
# import ABQ-packages
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
from textRepr import *

# import python packages
import numpy as np 
import sys
import os

################################################################################
# ASSIGN INPUT PARAMETERS AND DEFAULT PARAMETERS
################################################################################
# Default parameters
Volumenvernetzung               = True     # Volumenvernetzung generell ein/aus              | Default: True
Innenvernetzung                 = True     # Innenvernetzung generell ein/aus                | Default: True
Innenvernetzung_vorher          = False    # Zuerst Innenvernetzung ansonsten danach         | Default: False
Seed_gegenueberliegende_Seite   = False    # Seed x1 Seite                                   | Default: False
on_or_off                       = OFF      # Projektion bei Koordinantenkorrektur            | Default: OFF
Partnersets_erzeugen            = False    # erzeugen von Sets der Knoten ohne Partner       | Default: False

# definition of numerical toleranz for BoundingBox and validation function
DELTA=10.0**-7              # tolerance for BoundingBox
DELTA_Partner=10.0**-12     # accuracy, for checking copied nodes

# parameters for seed density

input_arguments = sys.argv[-1]
# transform input_arguments into dictionary
input_args = {}

for arg in input_arguments.split():
    try:
        input_args.update( { arg.split('=')[0]:float(arg.split('=')[-1]) } )
    except:
        input_args.update( { arg.split('=')[0]:arg.split('=')[-1] } )

################################################################################
# MESH PARAMETERS
mesh_density_matrix = 0.09   # Default: 0.09
if 'mesh_density_matrix' in input_args.keys():
    mesh_density_matrix       = input_args['mesh_density_matrix']

mesh_density_yarn_ellipse = 0.05   # Default: 0.05
if 'mesh_density_yarn_ellipse' in input_args.keys():
    mesh_density_yarn_ellipse       = input_args['mesh_density_yarn_ellipse']

mesh_density_yarn_inside = 0.06   # Default: 0.06
if 'mesh_density_yarn_inside' in input_args.keys():
    mesh_density_yarn_inside       = input_args['mesh_density_yarn_inside']

cut = False
Ausschnittsshift = True
Volumenvernetzung = True
# Modell erzeugen
mod_nam= 'RVE'
RVE = mdb.Model(name=mod_nam) 
if 'Model-1' in mdb.models:
    del mdb.models['Model-1']

################################################################################
# geometry parameters

# geometry yarn cross section
h  = input_args['h']  # height of ellipse
b  = input_args['b']  # width of ellipse  
d1 = input_args['d1'] # yarn distance in y-direction at maximum
T  = input_args['T']  # period length 
e  = 2*d1+2*h	      # thickness of matrix
d2 = T/2.0            # yarn distance in z-direction/x-direction

# Geometrie Matrix                        
x = T                # x-Richtung                        
y = e                # y-Richtung                        
z = T                # z-Richtung                        

# geometry of yarn-profile: 
splinesupports = 700   # number of splinesupports
A = d1/2.0 + h/2.0     # amplitude                               
n = 2                  # number of periods, n>1

sinuspoints = []    # list of coordinates of splinesupports
for i in range(splinesupports+1):                              # +1 to get odd number of splinesupports
   x_sin_coo = (float(i)/splinesupports)*T*n - (n*T)/2.0       # coordinate system located in middle of yarn
   y_sin_coo = (float(A))*sin(n*2*pi*(float(i)/splinesupports))
   sinuspoints.append([x_sin_coo, y_sin_coo])

################################################################################
# generate part of yarn and matrix
################################################################################
   
# Part: generate yarn
Pfad_sketch = RVE.ConstrainedSketch(name='Pfad', sheetSize=4.0) 
Pfad_sketch.Spline(points=sinuspoints) #Definiert Spline
Querschnitt_sketch = RVE.ConstrainedSketch(name='Querschnitt', sheetSize=4.0, transform=(-0.293692389017689, -0.955899984639127, 0.0, 0.0, 0.0, 1.0, -0.955899984639127, 0.293692389017689, 0.0, -4.0, 1.5, 0.0)) # Dreht Koordinatensystem sodass Querschnitt des Splines gezeichnet werden kann
Querschnitt_sketch.EllipseByCenterPerimeter(axisPoint1=(0.0, b/2.0), axisPoint2=(-h/2.0, 0.0), center=(0.0, 0.0)) # Ellipsenquerschnitt definieren

Faser_part = RVE.Part(dimensionality=THREE_D, name='Faser', type=DEFORMABLE_BODY)
Faser_part.BaseSolidSweep(path= RVE.sketches['Pfad'], sketch=RVE.sketches['Querschnitt'])
del RVE.sketches['Querschnitt']
del RVE.sketches['Pfad']   
faser = RVE.parts['Faser'] 

# Part: Matrix erzeugen
Matrix_sketch = RVE.ConstrainedSketch(name='Matrix_sketch', sheetSize=2.0)
Matrix_sketch.rectangle(point1=(0.0, 0.0),point2=(1.0*x, y))
Matrix_part = RVE.Part(dimensionality=THREE_D, name='Matrix', type=DEFORMABLE_BODY)
Matrix_part.BaseSolidExtrude(depth=1.0*z, sketch=Matrix_sketch)
del RVE.sketches['Matrix_sketch']
matrix = RVE.parts['Matrix'] 

#########Zuschneiden der Faser auf eine Periodendauer##########################
######### sorgt fuer parallel zu Koordinatensystemebenen liegenden Abschnnitt##

# Koordinatensytem Achsen und Ebenen
x_Achse  = faser.DatumAxisByPrincipalAxis(principalAxis=XAXIS)
y_Achse  = faser.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
z_Achse  = faser.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
xy_Ebene = faser.DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XYPLANE)
yz_Ebene = faser.DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=YZPLANE)
xz_Ebene = faser.DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XZPLANE)

# Bezugsebenen fuer den Schnitt in Nullstelle
Schnittebene_2 = faser.DatumPlaneByPrincipalPlane(offset=(1.0/2.0)*T, principalPlane=YZPLANE)
Schnittebene_1 = faser.DatumPlaneByPrincipalPlane(offset=-(1.0/2.0)*T, principalPlane=YZPLANE)

# Bezugsachsen fuer Schnitt (faser.edges[2] referenziert die aeusserste Ellipsenkante)
Edge1 = faser.DatumAxisByNormalToPlane(plane=faser.datums[xz_Ebene.id], point=faser.InterestingPoint(faser.edges[2], CENTER))
Edge2 = faser.DatumAxisByNormalToPlane(plane=faser.datums[xz_Ebene.id], point=faser.InterestingPoint(faser.edges[0], CENTER))

# Faserschnitt 1
RVE.ConstrainedSketch(gridSpacing=0.33, name='__profile__', 
    sheetSize=13.21, transform=
    faser.MakeSketchTransform(
    sketchPlane=faser.datums[Schnittebene_1.id], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=faser.datums[Edge1.id], 
    sketchOrientation=RIGHT, origin=(-0.0, 0.0, -0.0)))
faser.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=RVE.sketches['__profile__'])
# Rechteck fuer Cut
RVE.sketches['__profile__'].rectangle(point1=(-4.0*b, 2.0*A), point2=(4.0*b, -2.0 *A))
faser.CutExtrude(flipExtrudeDirection=OFF, sketch=
    RVE.sketches['__profile__'], sketchOrientation=RIGHT, 
    sketchPlane=faser.datums[Schnittebene_1.id], sketchPlaneSide=
    SIDE1, sketchUpEdge=faser.datums[Edge1.id])
del RVE.sketches['__profile__']

# Faserschnitt 2
RVE.ConstrainedSketch(gridSpacing=0.33, name='__profile__', 
    sheetSize=13.21, transform=
    faser.MakeSketchTransform(
    sketchPlane=faser.datums[Schnittebene_2.id], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=faser.datums[Edge2.id], 
    sketchOrientation=RIGHT, origin=(0.0, 0.0, -0.0)))
faser.projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=RVE.sketches['__profile__'])
# Rechteck fuer Cut
RVE.sketches['__profile__'].rectangle(point1=(-4.0*b, 2.0*A), point2=(4.0*b, -2.0*A))
faser.CutExtrude(flipExtrudeDirection=ON, sketch=
    RVE.sketches['__profile__'], sketchOrientation=RIGHT, 
    sketchPlane=faser.datums[Schnittebene_2.id], sketchPlaneSide=
    SIDE1, sketchUpEdge=faser.datums[Edge2.id])
del RVE.sketches['__profile__']

# Horizontale Partition der Faser erstellen
faser.PartitionCellByExtrudeEdge(cells= faser.cells[0], edges=(faser.edges[2], ),line=faser.datums[z_Achse.id], sense=REVERSE)

# Horizontale Partition der Faser erstellen
# faser.PartitionCellByExtrudeEdge(cells= faser.cells[0], edges=(faser.edges[1], ),line=faser.datums[z_Achse.id], sense=REVERSE)

####################################################################################################################
#####################################Zusammenbau der Einheitszelle##################################################
####################################################################################################################
Asm = RVE.rootAssembly

# Erstellen von Instanz fuer Matrix und Fasern

List_of_Instances = [
Asm.Instance(dependent=ON, name='KettFaser-1', part= faser),
Asm.Instance(dependent=ON, name='KettFaser-2', part= faser),
Asm.Instance(dependent=ON, name='SchussFaser-1', part= faser),
Asm.Instance(dependent=ON, name='SchussFaser-2', part= faser),
]
Asm.Instance(dependent=ON, name='Matrix-1', part= matrix)
# Anpassen der Koordinatensysteme durch Verschieben aller Fasern 
#um T/2 in x-Richtung und y/2 in y-Richtung sowie in z-Richtung um den halben Faserabstand
Asm.translate(instanceList=('KettFaser-1','KettFaser-2','KettFaser-3','KettFaser-4','KettFaser-5','KettFaser-6','SchussFaser-1','SchussFaser-2','SchussFaser-3','SchussFaser-4','SchussFaser-5','SchussFaser-6', ), vector=(T/2.0,y/2.0, d2/2.0))
Asm.translate(instanceList=('Matrix-1', ), vector=((-3.0/4.0)*T,0,(-3.0/4.0)*T) )

#Rotieren von Kettfaser um 180 Grad
Asm.rotate(angle=180.0, axisDirection=(1.0, 0.0, 0.0), axisPoint=(T/2.0,y/2.0, T/2.0), instanceList=('KettFaser-2','KettFaser-4','KettFaser-6',))
# Rotieren aller Schussfasern um 90 Grad
Asm.rotate(angle=90.0, axisDirection=(0.0, 1.0, 0.0), axisPoint=(T/2.0,y/2.0, T/2.0), instanceList=('SchussFaser-1', 'SchussFaser-2','SchussFaser-3', 'SchussFaser-4','SchussFaser-5','SchussFaser-6', ))
#Rotieren von Schussfaser um 180 Grad
Asm.rotate(angle=180.0, axisDirection=(0.0, 0.0, 1.0), axisPoint=(T/2.0,y/2.0, T/2.0), instanceList=('SchussFaser-2','SchussFaser-4','SchussFaser-6', ))


#neue Fasern verschieben
Asm.translate(instanceList=('KettFaser-3',), vector=(0.0,0.0, T))
Asm.translate(instanceList=('KettFaser-4',), vector=(0.0,0.0, -T))
Asm.translate(instanceList=('KettFaser-5',), vector=(0.0,0.0, -T))
Asm.translate(instanceList=('KettFaser-6',), vector=(0.0,0.0, T))
Asm.translate(instanceList=('SchussFaser-3',), vector=(T,0.0, 0.0))
Asm.translate(instanceList=('SchussFaser-4',), vector=(-T,0.0, 0.0))
Asm.translate(instanceList=('SchussFaser-5',), vector=(-T,0.0, 0.0))
Asm.translate(instanceList=('SchussFaser-6',), vector=(T,0.0, 0.0))

# koosys anpassen 
Asm.translate(instanceList=('SchussFaser-1', 'SchussFaser-2','SchussFaser-3', 'SchussFaser-4','SchussFaser-5','SchussFaser-6','KettFaser-1', 'KettFaser-2','KettFaser-3', 'KettFaser-4','KettFaser-5','KettFaser-6', 'Matrix-1', ), vector=(-(3.0/4.0)*T ,-0.5*y,-(3.0/4.0)*T))

#Matrix auf Fasern
Asm.translate(instanceList=('Matrix-1', ), vector=(T ,0,T))
# Koosy anpassen
Asm.translate(instanceList=('SchussFaser-1', 'SchussFaser-2','KettFaser-1', 'KettFaser-2','Matrix-1', ), vector=((1.0/2.0)*T ,0.5*y,(1.0/2.0)*T))


############MATERIAL##############
RVE.Material(name='Yarn-Material')
RVE.materials['Yarn-Material'].Elastic(table=((4350.0, 0.36),  ))
RVE.HomogeneousSolidSection(material='Yarn-Material',   name= 'Yarn-Section',  thickness=None)

if Ausschnittsshift == True:
    for Element in List_of_Instances:
        Asm.translate(instanceList = (Element.name, ), vector=(T/4.0,0.0, T/4.0))
                

if cut == True:

# Defining Cuttingplane, offset is any component of any translationvector 
    Cuttingplane_1a = unitcell.DatumPlaneByPrincipalPlane(offset= 0, principalPlane=XYPLANE)
    Cuttingplane_1b = unitcell.DatumPlaneByPrincipalPlane(offset= T, principalPlane=XYPLANE)
    Cuttingplane_2a = unitcell.DatumPlaneByPrincipalPlane(offset= 0, principalPlane=YZPLANE)
    Cuttingplane_2b = unitcell.DatumPlaneByPrincipalPlane(offset= T, principalPlane=YZPLANE)
    Cuttingplane_3a = unitcell.DatumPlaneByPrincipalPlane(offset= 0, principalPlane=XZPLANE)
    Cuttingplane_3b = unitcell.DatumPlaneByPrincipalPlane(offset= e, principalPlane=XZPLANE)
    # Defining ReferenceAxis needed for Cut
    ReferenceAxis1 = unitcell.DatumAxisByPrincipalAxis(principalAxis=XAXIS)
    ReferenceAxis2 = unitcell.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
    ReferenceAxis3 = unitcell.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)

    Cut = [(Cuttingplane_1a, ReferenceAxis1, OFF),
           (Cuttingplane_1b, ReferenceAxis1, ON),
           (Cuttingplane_2a, ReferenceAxis2, OFF),
           (Cuttingplane_2b, ReferenceAxis2, ON),
           (Cuttingplane_3a, ReferenceAxis3, OFF),
           (Cuttingplane_3b, ReferenceAxis3, ON)]

    # #### Create Cut########
    for Cuttingplane, Axis, Flip in Cut:
        RVE.ConstrainedSketch(gridSpacing=6.64, name='__profile__', sheetSize=265.6, transform=
            unitcell.MakeSketchTransform(
            sketchPlane=unitcell.datums[Cuttingplane.id],
            sketchPlaneSide=SIDE1, 
            sketchUpEdge=unitcell.datums[Axis.id], sketchOrientation=RIGHT, origin=(-0.0, 0.0, -0.0)))
        unitcell.projectReferencesOntoSketch(filter=
            COPLANAR_EDGES, sketch=RVE.sketches['__profile__'])
        RVE.sketches['__profile__'].rectangle(point1=(4.0*T, -4.0*T), point2=(-4.0*T, 4.0*T))
        unitcell.CutExtrude(flipExtrudeDirection=Flip, sketch=
            RVE.sketches['__profile__'], sketchOrientation=RIGHT, 
            sketchPlane=unitcell.datums[Cuttingplane.id], sketchPlaneSide=
            SIDE1, sketchUpEdge=unitcell.datums[Axis.id])
        del RVE.sketches['__profile__']

##################################################PERIODISCHES NETZ#########################################################

# Ausschneiden der Matrix
Asm.InstanceFromBooleanCut(cuttingInstances=(Asm.instances['KettFaser-1'],
                                             Asm.instances['KettFaser-2'],
                                             Asm.instances['SchussFaser-1'],
                                             Asm.instances['SchussFaser-2'], ), 
                           instanceToBeCut=Asm.instances['Matrix-1'], 
                           name='Matrix_with_Cut', 
                           originalInstances=SUPPRESS)

# vorherige Instanzen wieder Einfuegen
Asm.resumeFeatures(('KettFaser-1','KettFaser-2','SchussFaser-1', 'SchussFaser-2'))
                    
# Zusammenfuegen Instanzen zu einem Part
Asm.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(Asm.instances['KettFaser-1'], 
               Asm.instances['KettFaser-2'], 
               Asm.instances['SchussFaser-1'], 
               Asm.instances['SchussFaser-2'], 
               Asm.instances['Matrix_with_Cut-1']), 
    keepIntersections=ON,
    name='Unitcell',
    originalInstances=SUPPRESS)    
Asm.regenerate()                             
                             
unitcell=RVE.parts['Unitcell']
                          
# Koordinantesystem Anzeigen durch Bezugsachsen
unitcell.DatumAxisByPrincipalAxis(principalAxis=XAXIS)
unitcell.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
unitcell.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)



##########################################################################################################################################

####################################################################################################################
########################################## Mesh Control Einstellungen ##############################################
####################################################################################################################
#Lineares Netz
unitcell.setMeshControls(elemShape=TET, regions=unitcell.cells[:], technique=FREE)
unitcell.setElementType(elemTypes=(ElemType(elemCode=C3D8R, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
    elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)), 
    regions=(unitcell.cells[:], ))
    
####################################################################################################################
#############################Definition Boundingboxen fuer Aussenflaechen Matrix im Partmodus########################
####################################################################################################################
#BoundingBoxen
# Alle Kanten auf Ebene x0 mit BoundingBox in Alle_Kanten_x0 speichern
# Alle Kanten auf Ebene x1 mit BoundingBox in Alle_Kanten_x1 speichern
# ...
Alle_Kanten_x0 = unitcell.edges.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                 xMax=+DELTA,yMax=y+DELTA,zMax=z+DELTA)                                           
Alle_Kanten_x1 = unitcell.edges.getByBoundingBox(xMin=-DELTA+x,yMin=-DELTA,zMin=-DELTA,
                                                 xMax=+DELTA+x,yMax=y+DELTA,zMax=z+DELTA)                                              
Alle_Kanten_y0 = unitcell.edges.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                 xMax=+DELTA+x,yMax=+DELTA,zMax=+DELTA+z)                                           
Alle_Kanten_y1 = unitcell.edges.getByBoundingBox(xMin=-DELTA,yMin=-DELTA+y,zMin=-DELTA,
                                                 xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)
Alle_Kanten_z0 = unitcell.edges.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                 xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA)                                           
Alle_Kanten_z1 = unitcell.edges.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA+z,
                                                 xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)
# Alle Flaechen auf Ebene x0 mit BoundingBox in Alle_Flaechen_x0 speichern
# Alle Flaechen auf Ebene x1 mit BoundingBox in Alle_Flaechen_x1 speichern
# ...
Alle_Flaechen_x0 = unitcell.faces.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                 xMax=+DELTA,yMax=y+DELTA,zMax=z+DELTA)                                           
Alle_Flaechen_x1 = unitcell.faces.getByBoundingBox(xMin=-DELTA+x,yMin=-DELTA,zMin=-DELTA,
                                                 xMax=+DELTA+x,yMax=y+DELTA,zMax=z+DELTA)
Alle_Flaechen_y0 = unitcell.faces.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                 xMax=+DELTA+x,yMax=+DELTA,zMax=+DELTA+z)                                           
Alle_Flaechen_y1 = unitcell.faces.getByBoundingBox(xMin=-DELTA,yMin=-DELTA+y,zMin=-DELTA,
                                                 xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)
Alle_Flaechen_z0 = unitcell.faces.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                 xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA)                                           
Alle_Flaechen_z1 = unitcell.faces.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA+z,
                                                 xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)
                                                 
####################################################################################################################
########################################### Skript zum Zuordnen der Kanten zu den zugehoerigen Zellen ##############
####################################################################################################################
                                                                               
# Alle Kanten der jeweiligen Zelle im Dictionary Dic_Kanten = {}
# Alle Flaechen der jeweiligen Zelle im Dictionary Dic_Flaechen= {} speichern
Dic_Kanten   = {}          # Dictionary von allen Kanten, zugriff je Zelle ueber z.B. Dic_Kanten['Alle_Zellkanten1']
Dic_Flaechen = {}          # Dictionary von allen Flaechen, zugriff je Zelle ueber z.B. Dic_Flaechen['Zellflaeche1']
    
Dic_Kanten_in_x0    = {}   # Alle Kanten von Zelle in der x0-Ebene
Dic_Kanten_in_x1    = {}   # Alle Kanten von Zelle in der x1-Ebene
Dic_Kanten_in_y0    = {}   # Alle Kanten von Zelle in der y0-Ebene
Dic_Kanten_in_y1    = {}   # Alle Kanten von Zelle in der y1-Ebene
Dic_Kanten_in_z0    = {}   # Alle Kanten von Zelle in der z0-Ebene
Dic_Kanten_in_z1    = {}   # Alle Kanten von Zelle in der z1-Ebene
Dic_Kanten_inside   = {}   # Alle Kanten von Zelle innerhalb der Matrix

Dic_Flaechen_in_x0  = {}   # Alle Flaechen von Zelle in der x0-Ebene
Dic_Flaechen_in_x1  = {}   # Alle Flaechen von Zelle in der x1-Ebene
Dic_Flaechen_in_y0  = {}   # Alle Flaechen von Zelle in der y0-Ebene
Dic_Flaechen_in_y1  = {}   # Alle Flaechen von Zelle in der y1-Ebene
Dic_Flaechen_in_z0  = {}   # Alle Flaechen von Zelle in der z0-Ebene
Dic_Flaechen_in_z1  = {}   # Alle Flaechen von Zelle in der z1-Ebene
Dic_Flaechen_inside = {}   # Alle Flaechen von Zelle innerhalb der Matrix

# Fuer jede Zelle die Kanten finden und im Dictionary speichern
for i in range(len(unitcell.cells)):
    Zellkanten = []                                  # temporaere Liste mit Zellkanten, die fuer jeden Schleifendurchlauf ueberschrieben wird
    Kantenindizes = unitcell.cells[i].getEdges()     # Vektor mit Indexen der Kanten von Zelle i, z.B: (1, 7, 16, 17, 18, 19)
    for Index in Kantenindizes:
        Auswahlkante = unitcell.edges[Index]         # Vektor mit zugehoerigem Pfad zum Index, z.B. mdb.models['RVE'].parts['Unitcell'].edges.findAt((0.0, 0.155645, 1.47937))
        Zellkanten.append(Auswahlkante)
    Dic_Kanten.update({str(i):Zellkanten})           # Dictionary in dem die Eintraege von Zellkanten gespeichert weden

                                                 
    # Schnittmenge von Alle_Kanten_x0(gefunden mit BoundingBox) und allen Kanten der jeweiligen Zelle  bilden (Dic_Kanten[str(i)]) und Kanten einsortieren
    
    Kanten_Zelle_in_x0  = []   # Alle Kanten von Zelle in der x0-Ebene
    Kanten_Zelle_in_x1  = []   # Alle Kanten von Zelle in der x1-Ebene
    Kanten_Zelle_in_y0  = []   # Alle Kanten von Zelle in der y0-Ebene
    Kanten_Zelle_in_y1  = []   # Alle Kanten von Zelle in der y1-Ebene
    Kanten_Zelle_in_z0  = []   # Alle Kanten von Zelle in der z0-Ebene
    Kanten_Zelle_in_z1  = []   # Alle Kanten von Zelle in der z1-Ebene
    Kanten_Zelle_inside = []   # Alle Kanten von Zelle innerhalb der Matrix
    
    # Reihenfolge bei Sortierung wichtig, weil Ebenen ueberlappen
    for Kante in Dic_Kanten[str(i)]:
        if Kante in Alle_Kanten_x0:
            Kanten_Zelle_in_x0.append(Kante)
        elif Kante in Alle_Kanten_z0:
            Kanten_Zelle_in_z0.append(Kante)  
        elif Kante in Alle_Kanten_y0:
            Kanten_Zelle_in_y0.append(Kante)   
        elif Kante in Alle_Kanten_x1:
            Kanten_Zelle_in_x1.append(Kante)
        elif Kante in Alle_Kanten_y1:
            Kanten_Zelle_in_y1.append(Kante)
        elif Kante in Alle_Kanten_z1:
            Kanten_Zelle_in_z1.append(Kante) 
        else:
            Kanten_Zelle_inside.append(Kante)
            
    Dic_Kanten_in_x0.update({str(i):Kanten_Zelle_in_x0}) 
    Dic_Kanten_in_x1.update({str(i):Kanten_Zelle_in_x1}) 
    Dic_Kanten_in_y0.update({str(i):Kanten_Zelle_in_y0}) 
    Dic_Kanten_in_y1.update({str(i):Kanten_Zelle_in_y1}) 
    Dic_Kanten_in_z0.update({str(i):Kanten_Zelle_in_z0}) 
    Dic_Kanten_in_z1.update({str(i):Kanten_Zelle_in_z1}) 
    Dic_Kanten_inside.update({str(i):Kanten_Zelle_inside})
    
    ####################################################################################################################
    #######################################Skript zum Flaechen finden ############ #####################################
    ####################################################################################################################

    # Fuer jede Zelle die Flaeche finden und im Dictionary speichern

    Zellflaechen = []                                      # temporaere Liste mit Zellflaechen, die fuer jeden Schleifendurchlauf ueberschrieben wird
    Flaechenindizes = unitcell.cells[i].getFaces()         # Vektor mit Indexen der Flaechen von Zelle i, z.B: (6, 7, 8, 11)
    for Index in Flaechenindizes:
        Auswahlflaeche = unitcell.faces[Index]             # Vektor mit zugehoerigem Pfad zum Index,z.B.: mdb.models['RVE'].parts['Unitcell'].faces.findAt((1.25, 0.2, 1.25),)
        Zellflaechen.append(Auswahlflaeche)
    Dic_Flaechen.update( {str(i):Zellflaechen} )           # Dictionary in dem die Eintraege von Zellflaechen gespeichert weden

               
    # Schnittmenge von Alle_Flaechen__x0(gefunden mit BoundingBox) und allen Flaechen der jeweiligen Zelle  bilden (Dic_Flaechen[str(i)]) und Flaechen einsortieren
    
    Flaechen_Zelle_in_x0  = []   # Alle Flaechen von Zelle in der x0-Ebene
    Flaechen_Zelle_in_x1  = []   # Alle Flaechen von Zelle in der x1-Ebene
    Flaechen_Zelle_in_y0  = []   # Alle Flaechen von Zelle in der y0-Ebene
    Flaechen_Zelle_in_y1  = []   # Alle Flaechen von Zelle in der y1-Ebene
    Flaechen_Zelle_in_z0  = []   # Alle Flaechen von Zelle in der z0-Ebene
    Flaechen_Zelle_in_z1  = []   # Alle Flaechen von Zelle in der z1-Ebene
    Flaechen_Zelle_inside = []   # Alle Flaechen von Zelle innerhalb der Matrix
    
    for Flaeche in Dic_Flaechen[str(i)]:
        if Flaeche in Alle_Flaechen_x0:
            Flaechen_Zelle_in_x0.append(Flaeche)  
        if Flaeche in Alle_Flaechen_y0:
            Flaechen_Zelle_in_y0.append(Flaeche)
        if Flaeche in Alle_Flaechen_z0:
            Flaechen_Zelle_in_z0.append(Flaeche)    
        if Flaeche in Alle_Flaechen_x1:
            Flaechen_Zelle_in_x1.append(Flaeche)
        if Flaeche in Alle_Flaechen_y1:
            Flaechen_Zelle_in_y1.append(Flaeche)
        if Flaeche in Alle_Flaechen_z1:
            Flaechen_Zelle_in_z1.append(Flaeche) 
        else:
            if not Flaeche in Alle_Flaechen_x0:
                if not Flaeche in Alle_Flaechen_y0:
                    if not Flaeche in Alle_Flaechen_z0:
                        if not Flaeche in Alle_Flaechen_x1:
                            if not Flaeche in Alle_Flaechen_y1:
                                if not Flaeche in Alle_Flaechen_z1:
                                    Flaechen_Zelle_inside.append(Flaeche)
    
    Dic_Flaechen_in_x0.update({str(i):Flaechen_Zelle_in_x0}) 
    Dic_Flaechen_in_x1.update({str(i):Flaechen_Zelle_in_x1}) 
    Dic_Flaechen_in_y0.update({str(i):Flaechen_Zelle_in_y0}) 
    Dic_Flaechen_in_y1.update({str(i):Flaechen_Zelle_in_y1}) 
    Dic_Flaechen_in_z0.update({str(i):Flaechen_Zelle_in_z0}) 
    Dic_Flaechen_in_z1.update({str(i):Flaechen_Zelle_in_z1}) 
    Dic_Flaechen_inside.update({str(i):Flaechen_Zelle_inside})

    ####################################################################################################################
    ######################## Seed der Kanten ###########################################################################
    ####################################################################################################################
    #Seed der Matrix
    if len(unitcell.cells[i].getFaces()) > 4:                     # 4 ist die Anzahl der Flaechen einer Faser, Matrix 
                                                                  # hat mehr Flaechen und wird darueber zugeordnet
        Matrix_Aussenkanten = (unitcell.edges.findAt((x/2.0,0.0,0.0)),
                               unitcell.edges.findAt((x,y/2.0,0.0)),
                               unitcell.edges.findAt((x/2.0,y,0.0)),
                               unitcell.edges.findAt((0.0,y/2,0.0)),
                               unitcell.edges.findAt((0.0,0.0,z/2.0)),
                               unitcell.edges.findAt((0.0,y,z/2.0)),
                               unitcell.edges.findAt((0.0,y/2.0,z)))
        unitcell.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1, 
                            edges= (Matrix_Aussenkanten), size=mesh_density_matrix)
    else: 
        # Seed der Zelle innerhalb der Matrix
        unitcell.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1, 
                                edges= Kanten_Zelle_inside, size=mesh_density_yarn_inside)
        # Seed der Zelle Seite x0,y0,z0
        unitcell.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1, 
                                edges= Kanten_Zelle_in_x0, size=mesh_density_yarn_ellipse)
        unitcell.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1, 
                                edges= Kanten_Zelle_in_y0, size=mesh_density_yarn_ellipse)
        unitcell.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1, 
                                edges= Kanten_Zelle_in_z0, size=mesh_density_yarn_ellipse)
        if Seed_gegenueberliegende_Seite ==  True:
            #Seed der Zelle Seite x1 y1 und z1
            unitcell.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1, 
                                    edges= Kanten_Zelle_in_x1, size=mesh_density_yarn_ellipse)
            unitcell.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1, 
                                    edges= Kanten_Zelle_in_y1, size=mesh_density_yarn_ellipse)
            unitcell.seedEdgeBySize(constraint=FIXED, deviationFactor=0.01, 
                                    edges= Kanten_Zelle_in_z1, size=mesh_density_yarn_ellipse)

                 
####################################################################################################################
####################Vernetzen der Oberflaeche, Offset des Meshs, Assoziieren der Geometrie #########################
####################################################################################################################
# nicht einzeln moeglich, da auf Offsetelemente nur indirekt ueber Ausschlussverfahren zugegriffen werden kann
schon_zugeordnet = []               # Liste, welche schon zugeordnete Elementfaces enthaelt

if Innenvernetzung == True:
    if Innenvernetzung_vorher == True:
        # Innenliegende Faserflaechen vernetzen
        for i in range(len(unitcell.cells)):
            if not len(Dic_Flaechen_inside[str(i)]) == 0:
                unitcell.generateMesh(boundaryPreview=ON, regions= Dic_Flaechen_inside[str(i)], )
            else:
                pass
        schon_zugeordnet.extend(unitcell.elementFaces)  
            
# Fuer jede Zelle die Oberflaechen in x0,y0,z0-Ebenen vernetzen
for i in range(len(unitcell.cells)):
    print(i)
    # Vernetzen Oberflaeche von Zelle i in der x0-Ebene 
    if not len(Dic_Flaechen_in_x0[str(i)]) == 0:                                                                    # ohne diese Abfrage Fehlermeldung: list index out of range
        # Generiert Oberflachennetz der Zelle in x0-Ebene
        unitcell.generateMesh(boundaryPreview=ON, regions= Dic_Flaechen_in_x0[str(i)], )                            # Oberflaechenvernetzung von Zelle i in x0
        # Offset
        unitcell.generateMeshByOffset(distanceBetweenLayers=0.0, initialOffset=-x, meshType=SHELL, numLayers=1,     # Uebertragen des Netzes der Zelle i in x1-Ebene durch MeshByOffset
        offsetDirection=OUTWARD, region=Region(side1Elements=Dic_Flaechen_in_x0[str(i)][0].getElements()))          # Elemente sind die Meshelemente der x0 Flaeche der Zelle
        unitcell.mergeNodes(nodes= unitcell.nodes[:], tolerance=1e-8)                                               # alle x0-Elementfaces der Zelle, die aussortiert werden sollen
        ElementfacesX0=[]
        for j in range(len(Dic_Flaechen_in_x0[str(i)][0].getElements())):
            Elementflaechen = Dic_Flaechen_in_x0[str(i)][0].getElements()[j].getElemFaces()
            ElementfacesX0.extend(Elementflaechen)
        schon_zugeordnet.extend(ElementfacesX0) 
        #Speichern der zuzuweisenden Elemente ueber Ausschlussverfahren
        Elem_x1 = []                                                                                                # enthaelt alle Elementflaechen die zuzuweisen sind
        for Element in unitcell.elementFaces:
            if not Element in schon_zugeordnet:
                Elem_x1.append(Element)
        elem_x1 = MeshFaceArray(Elem_x1)                                                                        # elemFaces in associateMeshWithGeometry akzeptiert nur Datentyp MeshElementFaceArray oder Sets. Sets aus Elemfaces nicht erstellbar
        unitcell.associateMeshWithGeometry(elemFaces= elem_x1, geometricEntity=Dic_Flaechen_in_x1[str(i)][0])       # Geometriezuordnung
        schon_zugeordnet.extend(Elem_x1)                                                                            # alle jetzt zugeordneten Elemente aussortieren
    # # Vernetzen Oberflaeche von Zelle i in der y0-Ebene 
    if not len(Dic_Flaechen_in_y0[str(i)]) == 0:
        pass
        unitcell.generateMesh(boundaryPreview=ON, regions= Dic_Flaechen_in_y0[str(i)], )                            # Oberflaechenvernetzung von Zelle i in y0
        unitcell.generateMeshByOffset(distanceBetweenLayers=0.0, initialOffset=-y, meshType=SHELL, numLayers=1,     # Uebertragen des Netzes der Zelle i in x1-Ebene durch MeshByOffset 
        offsetDirection=OUTWARD, region=Region(side1Elements=Dic_Flaechen_in_y0[str(i)][0].getElements()))          # Elemente sind die Meshelemente der x0 Flaeche der Zelle
        unitcell.mergeNodes(nodes= unitcell.nodes[:], tolerance=1e-8)
        ElementfacesY0=[]
        for j in range(len(Dic_Flaechen_in_y0[str(i)][0].getElements())):
            Elementflaechen = Dic_Flaechen_in_y0[str(i)][0].getElements()[j].getElemFaces()
            ElementfacesY0.extend(Elementflaechen)
        schon_zugeordnet.extend(ElementfacesY0)                                                                     # alle y0-Elementfaces der Zelle, die aussortiert werden sollen
        Elem_y1 = []
        for Element in unitcell.elementFaces:
            if not Element in schon_zugeordnet:
                Elem_y1.append(Element)
        elem_y1 = MeshFaceArray(Elem_y1)                                                                        # elemFaces in associateMeshWithGeometry akzeptiert nur Datentyp MeshElementFaceArray oder Sets. Sets aus Elemfaces nicht erstellbar
        unitcell.associateMeshWithGeometry(elemFaces= elem_y1, geometricEntity=Dic_Flaechen_in_y1[str(i)][0])       # Geometriezuordnung
        schon_zugeordnet.extend(Elem_y1)          
    # Vernetzen Oberflaeche von Zelle i in der z0-Ebene 
    if not len(Dic_Flaechen_in_z0[str(i)]) == 0:
        unitcell.generateMesh(boundaryPreview=ON, regions= Dic_Flaechen_in_z0[str(i)], )                            # Oberflaechenvernetzung von Zelle i in z0
        unitcell.generateMeshByOffset(distanceBetweenLayers=0.0, initialOffset=-z, meshType=SHELL, numLayers=1,     # Uebertragen des Netzes der Zelle i in z1-Ebene durch MeshByOffset
        offsetDirection=OUTWARD, region=Region(side1Elements=Dic_Flaechen_in_z0[str(i)][0].getElements()))          # Elemente sind die Meshelemente der z0 Flaeche der Zelle
        unitcell.mergeNodes(nodes= unitcell.nodes[:], tolerance=1e-8)
        ElementfacesZ0=[]
        for j in range(len(Dic_Flaechen_in_z0[str(i)][0].getElements())):
            Elementflaechen = Dic_Flaechen_in_z0[str(i)][0].getElements()[j].getElemFaces()
            ElementfacesZ0.extend(Elementflaechen)                                   # alle z0-Elementfaces der Zelle, die aussortiert werden sollen
        schon_zugeordnet.extend(ElementfacesZ0)
        Elem_z1 = []
        for Element in unitcell.elementFaces:
            if not Element in schon_zugeordnet:
                Elem_z1.append(Element)
        elem_z1 = MeshFaceArray(Elem_z1)                                                                        # elemFaces in associateMeshWithGeometry akzeptiert nur Datentyp MeshElementFaceArray oder Sets. Sets aus Elemfaces nicht erstellbar
        unitcell.associateMeshWithGeometry(elemFaces= elem_z1, geometricEntity=Dic_Flaechen_in_z1[str(i)][0])       # Geometriezuordnung
        schon_zugeordnet.extend(Elem_z1)

####Speichern der Knotenkoordinaten vor dem Innenvernetzen, um zu vergleichen, ob diese sich innerhalb einer Ebene durch das Innenvernetzen veraendern         
# Koordinanten in x0 x1 vor dem Innenvernetzen
all_nodes_x0_before = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA, xMax=+DELTA,yMax=+DELTA+y,zMax=+DELTA+z )
all_nodes_x1_before = unitcell.nodes.getByBoundingBox(xMin=-DELTA+x,yMin=-DELTA,zMin=-DELTA,xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z )                                                                         
x0_coo_before = np.array([Knoten.coordinates for Knoten in all_nodes_x0_before])
x1_coo_before = np.array([Knoten.coordinates for Knoten in all_nodes_x1_before])
# Koordinanten in y0 y1 vor dem Innenvernetzen
all_nodes_y0_before = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,xMax=+DELTA+x,yMax=+DELTA,zMax=+DELTA+z)
all_nodes_y1_before = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA+y,zMin=-DELTA,xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)                                                                         
y0_coo_before = np.array([Knoten.coordinates for Knoten in all_nodes_y0_before])
y1_coo_before = np.array([Knoten.coordinates for Knoten in all_nodes_y1_before])
# Koordinanten in z0 z1 vor dem Innenvernetzen
all_nodes_z0_before = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA)
all_nodes_z1_before = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA+z,xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)                                                                        
z0_coo_before = np.array([Knoten.coordinates for Knoten in all_nodes_z0_before])
z1_coo_before = np.array([Knoten.coordinates for Knoten in all_nodes_z1_before])  

# Innenvernetzung: Vernetzung der Oberflachen der innerhalb der Matrix liegenden Faeden
if Innenvernetzung == True:
    if Innenvernetzung_vorher == False: 
        #Innenliegende Faserflaechen vernetzen
        for i in range(len(unitcell.cells)):
            unitcell.generateMesh(boundaryPreview=ON, regions= Dic_Flaechen_inside[str(i)], )      
    print('Innenvernetzung')

### Speichern der Knotenkoordinaten nach dem Innenvernetzen, um zu vergleichen, ob diese sich innerhalb einer Ebene durch das Innenvernezten veraendern    
# Koordinanten in x0 x1 nach dem Innenvernetzen
all_nodes_x0_after = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA, xMax=+DELTA,yMax=+DELTA+y,zMax=+DELTA+z )
all_nodes_x1_after = unitcell.nodes.getByBoundingBox(xMin=-DELTA+x,yMin=-DELTA,zMin=-DELTA,xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z )                                                                         
x0_coo_after = np.array([Knoten.coordinates for Knoten in all_nodes_x0_after])
x1_coo_after = np.array([Knoten.coordinates for Knoten in all_nodes_x1_after])
# Koordinanten in y0 y1 nach dem Innenvernetzen
all_nodes_y0_after = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,xMax=+DELTA+x,yMax=+DELTA,zMax=+DELTA+z)
all_nodes_y1_after = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA+y,zMin=-DELTA,xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)                                                                         
y0_coo_after = np.array([Knoten.coordinates for Knoten in all_nodes_y0_after])
y1_coo_after = np.array([Knoten.coordinates for Knoten in all_nodes_y1_after])
# Koordinanten in z0 z1 nach dem Innenvernetzen
all_nodes_z0_after = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA)
all_nodes_z1_after = unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA+z,xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)                                                                        
z0_coo_after = np.array([Knoten.coordinates for Knoten in all_nodes_z0_after])
z1_coo_after = np.array([Knoten.coordinates for Knoten in all_nodes_z1_after])  

## Ausgabe der Egebnisse des Koordinantenvergleichs: True: Alle Koordinanten sind gleich, False: mindestens ein Wert ist anders
print('Vergleich Koordinanten x0', np.all(x0_coo_before == x0_coo_after))    # np.all() prueft Wahrheit aller Werte des Arrays
print('Vergleich Koordinanten x1', np.all(x1_coo_before == x1_coo_after))
print('Vergleich Koordinanten y0', np.all(y0_coo_before == y0_coo_after))
print('Vergleich Koordinanten y1', np.all(y1_coo_before == y1_coo_after))
print('Vergleich Koordinanten z0', np.all(z0_coo_before == z0_coo_after))
print('Vergleich Koordinanten z1', np.all(z1_coo_before == z1_coo_after))

#Nachfolgend wird die Periodizitaet an verschiedenen Stelle ueberpruft: Hier: nach Oberflachennetz, dann nach Koordinantemanipulation
#Zuletzt nach Volumenvernetzung
####################################################################################################################
####################### Ueberpruefung der Periodischen Randbedinungen 1. Durchgang##################################
####################################################################################################################

#Set definieren mit Knoten von Flaeche X0 ueber BoundingBox, usw.
all_nodes_x0 = unitcell.Set(name='all_nodes_x0_1_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA,yMax=+DELTA+y,zMax=+DELTA+z ))
all_nodes_x1 = unitcell.Set(name='all_nodes_x1_1_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA+x,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z ))
all_nodes_y0 = unitcell.Set(name='all_nodes_y0_1_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA,zMax=+DELTA+z ))
all_nodes_y1 = unitcell.Set(name='all_nodes_y1_1_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA+y,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z))
all_nodes_z0 = unitcell.Set(name='all_nodes_z0_1_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA))
all_nodes_z1 = unitcell.Set(name='all_nodes_z1_1_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA+z,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z))                                                                                        

# Liste der Knotenkoordinaten von X0 mit list comprehension und numpy
x0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_x0.nodes])
y0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_y0.nodes])
z0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_z0.nodes])
x1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_x1.nodes])
y1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_y1.nodes])
z1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_z1.nodes])

#Vergleichsebene bestimmt zwei der Knotenkoordinaten fuer den Vergleich
YZ = [1,2] # fuer Ebene x
XZ = [0,2] # fuer Ebene y
XY = [0,1] # fuer Ebene z

#############Definition Validierungsfunktion###########################################
def finde_Knoten_ohne_Partner(Ausgangsebene_coo,Gegenueber_coo,Ebene_fuer_Abstandskontrolle,DELTA_Partner,all_nodes_Ausgangsebene):
    Kein_Partner = [] # Liste mit Knoten, die keinen Partner besitzen
    zaehler = 0  
    for i,Knotenkoordinate in enumerate(Ausgangsebene_coo):
        abstandsmatrix = Gegenueber_coo[:,Ebene_fuer_Abstandskontrolle] - Knotenkoordinate[Ebene_fuer_Abstandskontrolle]
        abstand1 = abstandsmatrix[:,0]
        abstand2 = abstandsmatrix[:,1]
        partner = min(np.sqrt(abstand1**2+abstand2**2)) < DELTA_Partner  # jeder Knoten hat einen Partner 
        if partner == False:
           zaehler = zaehler +1
           Kein_Partner.append(all_nodes_Ausgangsebene.nodes[i])
    return (Kein_Partner,zaehler)
# Rueckgabewerte: Kein_Partner: Liste  mit Knoten die keinen Partner besitzen
#                 zaheler: Wie viele Knoten keinen Partner besitzen
Kein_Partner_Seite_x0,zaehler_x  = finde_Knoten_ohne_Partner(x0_coo,x1_coo,YZ,DELTA_Partner,all_nodes_x0)    
Kein_Partner_Seite_y0,zaehler_y  = finde_Knoten_ohne_Partner(y0_coo,y1_coo,XZ,DELTA_Partner,all_nodes_y0)
Kein_Partner_Seite_z0,zaehler_z  = finde_Knoten_ohne_Partner(z0_coo,z1_coo,XY,DELTA_Partner,all_nodes_z0)
Kein_Partner_Seite_x1,zaehler_x1 = finde_Knoten_ohne_Partner(x1_coo,x0_coo,YZ,DELTA_Partner,all_nodes_x1)
Kein_Partner_Seite_y1,zaehler_y1 = finde_Knoten_ohne_Partner(y1_coo,y0_coo,XZ,DELTA_Partner,all_nodes_y1)
Kein_Partner_Seite_z1,zaehler_z1 = finde_Knoten_ohne_Partner(z1_coo,z0_coo,XY,DELTA_Partner,all_nodes_z1)

# Gibt laenge der Koordinatenlisten aus, um zu ueberpruefen, ob gleich viele Knoten auf 0 und 1 Ebene sind
print('1. Durchgang Anzahl der ueberprueften Koordinanten in x0 ', len(x0_coo))
print('1. Durchgang Anzahl der ueberprueften Koordinanten in x1 ', len(x1_coo))
print('1. Durchgang Anzahl der ueberprueften Koordinanten in y0 ', len(y0_coo))
print('1. Durchgang Anzahl der ueberprueften Koordinanten in y1 ', len(y1_coo))
print('1. Durchgang Anzahl der ueberprueften Koordinanten in z0 ', len(z0_coo))
print('1. Durchgang Anzahl der ueberprueften Koordinanten in z1 ', len(z1_coo))
# Gibt an, wie viele Knoten keinen Partner haben
print('1. Durchgang zaehler x',zaehler_x)
print('1. Durchgang zaehler y',zaehler_y)
print('1. Durchgang zaehler_z',zaehler_z)

if Partnersets_erzeugen == True:
    # Sets erzeugen um Punkte sichtbar zu machen
    if not len(Kein_Partner_Seite_x0)==0:
        unitcell.Set(name='kein_Partner_Seite_x0', nodes=MeshNodeArray(Kein_Partner_Seite_x0))
    if not len(Kein_Partner_Seite_x1)==0:    
        unitcell.Set(name='kein_Partner_Seite_x1', nodes=MeshNodeArray(Kein_Partner_Seite_x1))
    if not len(Kein_Partner_Seite_y0)==0:
        unitcell.Set(name='kein_Partner_Seite_y0', nodes=MeshNodeArray(Kein_Partner_Seite_y0))
    if not len(Kein_Partner_Seite_y1)==0:
        unitcell.Set(name='kein_Partner_Seite_y1', nodes=MeshNodeArray(Kein_Partner_Seite_y1))
    if not len(Kein_Partner_Seite_z0)==0:
        unitcell.Set(name='kein_Partner_Seite_z0', nodes=MeshNodeArray(Kein_Partner_Seite_z0))
    if not len(Kein_Partner_Seite_z1)==0:
        unitcell.Set(name='kein_Partner_Seite_z1', nodes=MeshNodeArray(Kein_Partner_Seite_z1))

#erzeugte Sets des vorherigen durchganges loeschen, da diese ansonsten als ungueltig (invalid) angezeigt werden
# Um Knoten ohne Partner des jeweiligen Durchgangs zu sehen: pause verwenden
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_x0_1_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_x1_1_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_y0_1_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_y1_1_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_z0_1_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_z1_1_Durchgang']
####################################################################################################################
####################### Koordinatenmanipulation#####################################################################
####################################################################################################################  
# Knoten der 1-Ebene werden Koordinaten der Knoten der 0-Ebene zugewiesen
# Projektion der Geometrie kann auf OFF gesetzt werden, da nur innerhalb der Ebene die Koordinaten verschoben werden
# Passende Partner finden ueber Minimum der Abstaende und Koordinaten anpassen (mit x0)
# Nur wenn fehlerhafte Knoten im ersten Durchgang

print('Koordinatenkorrektur, wenn Knoten ohne Partner im 1. Durchgang')
# Indizierung der Ebene, es werden nur die Abstaende von zwei der drei Koordinaten ueberpruft. 
# Die dritte Koordinate definiert den Abstand von 0 und 1 Ebene
# Gibt an, ob auf Spalte mit x-,y- oder z-Koordinaten des Arrays zugegriffen werden soll
X=0
Y=1
Z=2

# Fuer jeden Knoten aus der Liste der Knoten die keinen Partner haben wird auf der gegenueberliegenden Seite der Knoten mit 
# dem geringsten Abstand gewaehlt und die Koordinaten von 0 diesem zugewiesen
for Knoten_x0 in Kein_Partner_Seite_x0:                     
    Abstand=[]  # Liste der Abstaende aller Knoten aus x1 zu dem Bezugsknoten aus x0
    for j, Knoten_x1 in enumerate(Kein_Partner_Seite_x1):                           
        Knotenabstand = np.sqrt((Knoten_x0.coordinates[Y]-Knoten_x1.coordinates[Y])**2 +    # fuer jeden Knoten wird Abstand zu Knoten in x0 bestimmt
                          (Knoten_x0.coordinates[Z]-Knoten_x1.coordinates[Z])**2)
        Abstand.append((Knotenabstand, Knoten_x1,j))                                        # Abstand sowie zugehoeriger Knoten und durch j gekennzeichnete 
                                                                                            # Stelle im Vektor wird gespeichert
    Minimum = min(Abstand)                                                                  # Kleinster aller Abstaende wird in Minimum gespeichert
    Knoten = Minimum[1]                                                                     # zugehoriger Knoten wird gespeichert
    unitcell.editNode(coordinate2=Knoten_x0.coordinates[1],                                 # Koordinaten werden veraendert, Koordinaten von Knoten mit Minimalem Abstand zugewiesen
                      coordinate3=Knoten_x0.coordinates[2], 
                      nodes=Knoten, projectToGeometry=on_or_off)

for Knoten_y0 in Kein_Partner_Seite_y0:
    Abstand=[]
    for j, Knoten_y1 in enumerate(Kein_Partner_Seite_y1):
        Knotenabstand = np.sqrt((Knoten_y0.coordinates[X]-Knoten_y1.coordinates[X])**2 +
                          (Knoten_y0.coordinates[Z]-Knoten_y1.coordinates[Z])**2)
        Abstand.append((Knotenabstand, Knoten_y1,j))
    Minimum = min(Abstand)
    Knoten = Minimum[1]
    unitcell.editNode(coordinate1=Knoten_y0.coordinates[0], 
                      coordinate3=Knoten_y0.coordinates[2], 
                      nodes=Knoten, projectToGeometry=on_or_off)
                      
for Knoten_z0 in Kein_Partner_Seite_z0:
    Abstand=[]
    for j, Knoten_z1 in enumerate(Kein_Partner_Seite_z1):
        Knotenabstand = np.sqrt((Knoten_z0.coordinates[X]-Knoten_z1.coordinates[X])**2 +
                          (Knoten_z0.coordinates[Y]-Knoten_z1.coordinates[Y])**2)
        Abstand.append((Knotenabstand, Knoten_z1,j))
    Minimum = min(Abstand)
    Knoten = Minimum[1]
    unitcell.editNode(coordinate1=Knoten_z0.coordinates[0], 
                      coordinate2=Knoten_z0.coordinates[1], 
                      nodes=Knoten, projectToGeometry=on_or_off)


####################################################################################################################
####################### Ueberpruefung der Periodischen Randbedinungen 2. Durchgang##################################
####################################################################################################################


#Set definieren mit Knoten von Flaeche X0 ueber BoundingBox, usw.
all_nodes_x0 = unitcell.Set(name='all_nodes_x0_2_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA,yMax=+DELTA+y,zMax=+DELTA+z ))
all_nodes_x1 = unitcell.Set(name='all_nodes_x1_2_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA+x,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z ))
all_nodes_y0 = unitcell.Set(name='all_nodes_y0_2_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA,zMax=+DELTA+z ))
all_nodes_y1 = unitcell.Set(name='all_nodes_y1_2_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA+y,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z))
all_nodes_z0 = unitcell.Set(name='all_nodes_z0_2_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA))
all_nodes_z1 = unitcell.Set(name='all_nodes_z1_2_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA+z,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z))                                                                                        
                                                                                        
# Liste der Knotenkoordinaten von X0 mit list comprehension und numpy, usw.
x0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_x0.nodes])
y0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_y0.nodes])
z0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_z0.nodes])
x1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_x1.nodes])
y1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_y1.nodes])
z1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_z1.nodes])

Kein_Partner_Seite_x0,zaehler_x  = finde_Knoten_ohne_Partner(x0_coo,x1_coo,YZ,DELTA_Partner,all_nodes_x0)    
Kein_Partner_Seite_y0,zaehler_y  = finde_Knoten_ohne_Partner(y0_coo,y1_coo,XZ,DELTA_Partner,all_nodes_y0)
Kein_Partner_Seite_z0,zaehler_z  = finde_Knoten_ohne_Partner(z0_coo,z1_coo,XY,DELTA_Partner,all_nodes_z0)
Kein_Partner_Seite_x1,zaehler_x1 = finde_Knoten_ohne_Partner(x1_coo,x0_coo,YZ,DELTA_Partner,all_nodes_x1)
Kein_Partner_Seite_y1,zaehler_y1 = finde_Knoten_ohne_Partner(y1_coo,y0_coo,XZ,DELTA_Partner,all_nodes_y1)
Kein_Partner_Seite_z1,zaehler_z1 = finde_Knoten_ohne_Partner(z1_coo,z0_coo,XY,DELTA_Partner,all_nodes_z1)

print('2. Durchgang Anzahl der ueberprueften Koordinanten in x0 ', len(x0_coo))
print('2. Durchgang Anzahl der ueberprueften Koordinanten in x1 ', len(x1_coo))
print('2. Durchgang Anzahl der ueberprueften Koordinanten in y0 ', len(y0_coo))
print('2. Durchgang Anzahl der ueberprueften Koordinanten in y1 ', len(y1_coo))
print('2. Durchgang Anzahl der ueberprueften Koordinanten in z0 ', len(z0_coo))
print('2. Durchgang Anzahl der ueberprueften Koordinanten in z1 ', len(z1_coo))

print('2. Durchgang zaehler x',zaehler_x)
print('2. Durchgang zaehler y',zaehler_y)
print('2. Durchgang zaehler_z',zaehler_z)

if Partnersets_erzeugen == True:
# Sets erzeugen um Punkte sichtbar zu machen
    if not len(Kein_Partner_Seite_x0)==0:
        unitcell.Set(name='kein_Partner_Seite_x0_2_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_x0))
    if not len(Kein_Partner_Seite_x1)==0:    
        unitcell.Set(name='kein_Partner_Seite_x1_2_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_x1))
    if not len(Kein_Partner_Seite_y0)==0:
        unitcell.Set(name='kein_Partner_Seite_y0_2_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_y0))
    if not len(Kein_Partner_Seite_y1)==0:
        unitcell.Set(name='kein_Partner_Seite_y1_2_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_y1))
    if not len(Kein_Partner_Seite_z0)==0:
        unitcell.Set(name='kein_Partner_Seite_z0_2_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_z0))
    if not len(Kein_Partner_Seite_z1)==0:
        unitcell.Set(name='kein_Partner_Seite_z1_2_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_z1))

#erzeugte Sets des vorherigen durchganges loeschen
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_x0_2_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_x1_2_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_y0_2_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_y1_2_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_z0_2_Durchgang']
del mdb.models['RVE'].parts['Unitcell'].sets['all_nodes_z1_2_Durchgang']


####################################################################################################################
####################### Volumenvernetzung###########################################################################
####################################################################################################################

#Volumenvernetzung        
if Volumenvernetzung == True:
    print('Volumenvernetzung')
    unitcell.generateMesh()           # Volumenvernetzung

####################################################################################################################
####################### Ueberpruefung der Periodischen Randbedinungen 3. Durchgang##################################
####################################################################################################################


#Set definieren mit Knoten von Flaeche X0 ueber BoundingBox, usw.
all_nodes_x0 = unitcell.Set(name='all_nodes_x0_3_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA,yMax=+DELTA+y,zMax=+DELTA+z ))
all_nodes_x1 = unitcell.Set(name='all_nodes_x1_3_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA+x,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z ))
all_nodes_y0 = unitcell.Set(name='all_nodes_y0_3_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA,zMax=+DELTA+z ))
all_nodes_y1 = unitcell.Set(name='all_nodes_y1_3_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA+y,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z))
all_nodes_z0 = unitcell.Set(name='all_nodes_z0_3_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA))
all_nodes_z1 = unitcell.Set(name='all_nodes_z1_3_Durchgang', nodes= unitcell.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA+z,
                                                                                        xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z))                                                                                        
                                                                                        
# Liste der Knotenkoordinaten von X0 mit list comprehension und numpy, usw.
x0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_x0.nodes])
y0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_y0.nodes])
z0_coo = np.array([Knoten.coordinates for Knoten in all_nodes_z0.nodes])
x1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_x1.nodes])
y1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_y1.nodes])
z1_coo = np.array([Knoten.coordinates for Knoten in all_nodes_z1.nodes])

Kein_Partner_Seite_x0,zaehler_x  = finde_Knoten_ohne_Partner(x0_coo,x1_coo,YZ,DELTA_Partner,all_nodes_x0)    
Kein_Partner_Seite_y0,zaehler_y  = finde_Knoten_ohne_Partner(y0_coo,y1_coo,XZ,DELTA_Partner,all_nodes_y0)
Kein_Partner_Seite_z0,zaehler_z  = finde_Knoten_ohne_Partner(z0_coo,z1_coo,XY,DELTA_Partner,all_nodes_z0)
Kein_Partner_Seite_x1,zaehler_x1 = finde_Knoten_ohne_Partner(x1_coo,x0_coo,YZ,DELTA_Partner,all_nodes_x1)
Kein_Partner_Seite_y1,zaehler_y1 = finde_Knoten_ohne_Partner(y1_coo,y0_coo,XZ,DELTA_Partner,all_nodes_y1)
Kein_Partner_Seite_z1,zaehler_z1 = finde_Knoten_ohne_Partner(z1_coo,z0_coo,XY,DELTA_Partner,all_nodes_z1)

print('3. Durchgang Anzahl der ueberprueften Koordinanten in x0 ', len(x0_coo))
print('3. Durchgang Anzahl der ueberprueften Koordinanten in x1 ', len(x1_coo))
print('3. Durchgang Anzahl der ueberprueften Koordinanten in y0 ', len(y0_coo))
print('3. Durchgang Anzahl der ueberprueften Koordinanten in y1 ', len(y1_coo))
print('3. Durchgang Anzahl der ueberprueften Koordinanten in z0 ', len(z0_coo))
print('3. Durchgang Anzahl der ueberprueften Koordinanten in z1 ', len(z1_coo))

print('3. Durchgang zaehler x',zaehler_x)
print('3. Durchgang zaehler y',zaehler_y)
print('3. Durchgang zaehler_z',zaehler_z)

if Partnersets_erzeugen == True:
    # Sets erzeugen um Punkte sichtbar zu machen
    if not len(Kein_Partner_Seite_x0)==0:
        unitcell.Set(name='kein_Partner_Seite_x0_3_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_x0))
    if not len(Kein_Partner_Seite_x1)==0:    
        unitcell.Set(name='kein_Partner_Seite_x1_3_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_x1))
    if not len(Kein_Partner_Seite_y0)==0:
        unitcell.Set(name='kein_Partner_Seite_y0_3_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_y0))
    if not len(Kein_Partner_Seite_y1)==0:
        unitcell.Set(name='kein_Partner_Seite_y1_3_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_y1))
    if not len(Kein_Partner_Seite_z0)==0:
        unitcell.Set(name='kein_Partner_Seite_z0_3_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_z0))
    if not len(Kein_Partner_Seite_z1)==0:
        unitcell.Set(name='kein_Partner_Seite_z1_3_Durchgang', nodes=MeshNodeArray(Kein_Partner_Seite_z1))

        
    
####################################################################################################################
#############################ASSEMBLY MODUS#########################################################################
####################################################################################################################   
Asm.regenerate() # ansonsten sind Sets leer

Asm_inst = Asm.instances['Unitcell-1']

                                        
# #Set definieren mit Knoten von Flaeche X0 ueber BoundingBox, usw.
all_nodes_x0_Asm = Asm.Set(name='all_nodes_x0_Asm', nodes= Asm_inst.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                           xMax=+DELTA,yMax=+DELTA+y,zMax=+DELTA+z ))
all_nodes_x1_Asm = Asm.Set(name='all_nodes_x1_Asm', nodes= Asm_inst.nodes.getByBoundingBox(xMin=-DELTA+x,yMin=-DELTA,zMin=-DELTA,
                                                                                           xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z ))
all_nodes_y0_Asm = Asm.Set(name='all_nodes_y0_Asm', nodes= Asm_inst.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                           xMax=+DELTA+x,yMax=+DELTA,zMax=+DELTA+z ))
all_nodes_y1_Asm = Asm.Set(name='all_nodes_y1_Asm', nodes= Asm_inst.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA+y,zMin=-DELTA,
                                                                                           xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z))
all_nodes_z0_Asm = Asm.Set(name='all_nodes_z0_Asm', nodes= Asm_inst.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA,
                                                                                           xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA))
all_nodes_z1_Asm = Asm.Set(name='all_nodes_z1_Asm', nodes= Asm_inst.nodes.getByBoundingBox(xMin=-DELTA,yMin=-DELTA,zMin=-DELTA+z,
                                                                                           xMax=+DELTA+x,yMax=+DELTA+y,zMax=+DELTA+z)) 
####################################################################################################################
#######################Material ####################################################################################
####################################################################################################################
RVE.Material(name='Yarn-Material')
RVE.materials['Yarn-Material'].Elastic(table=((10590.0, 0.15),  ))

#### Matrix Material
RVE.Material(name='Matrix-Material')
RVE.materials['Matrix-Material'].Elastic(table=((4350.0, 0.36),  ))
####################################################################################################################
####################### Section definieren und Sets zuweisen #######################################################
####################################################################################################################

RVE.HomogeneousSolidSection(material='Yarn-Material',   name= 'Yarn-Section',  thickness=None)
RVE.HomogeneousSolidSection(material='Matrix-Material', name='Matrix-Section', thickness=None)

# Fasern
Zellen_pointON = []                     # fuer jede Zelle wird ein zugehoeriger Punkt definiert
for i in range(len(unitcell.cells)):
    if len(unitcell.cells[i].getFaces()) == 4:
        Zellen_pointON.append(unitcell.cells[i].pointOn)

unitcell.Set(cells= unitcell.cells.findAt(*Zellen_pointON), name='YARN') # *ermoeglicht an Argumente ranzukommen: Alle Punkte aus Zellen_pointON
unitcell.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region= unitcell.sets['YARN'], sectionName='Yarn-Section', thicknessAssignment=FROM_SECTION)

# Matrix
Zellen_pointON = []
for i in range(len(unitcell.cells)):
    if len(unitcell.cells[i].getFaces()) > 4:
        Zellen_pointON.append(unitcell.cells[i].pointOn)
        
unitcell.Set(cells= unitcell.cells.findAt(*Zellen_pointON), name='MATRIX')
unitcell.SectionAssignment(offset=0.0, offsetField= '', offsetType=MIDDLE_SURFACE, region=unitcell.sets['MATRIX'], sectionName='Matrix-Section', thicknessAssignment=FROM_SECTION)


# Job erzeugen
mdb.Job(activateLoadBalancing=False, atTime=None, contactPrint=OFF, 
    description='', echoPrint=OFF, explicitPrecision=SINGLE, 
    getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=
    PERCENTAGE, model='RVE', modelPrint=OFF, multiprocessingMode=DEFAULT, 
    name=input_args['inp_file_name'], 
    nodalOutputPrecision=SINGLE, numCpus=1, numDomains=1, 
    parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=ODB, 
    scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
# Inputfile erzeugen
mdb.jobs[input_args['inp_file_name']].writeInput()

# CLEAN ABQ-TEMP-FILES
for file in os.listdir(os.getcwd()):
    if file.startswith('abaqus.'):
        os.remove(file)
