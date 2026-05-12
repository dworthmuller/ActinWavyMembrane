
from dolfin import *
import numpy as np
import random
import copy
# import matplotlib.pyplot as plt
from ufl import tanh, cos
from numpy import savetxt 
import os
# import meshio
import gmsh
import pygmsh
import pandas as pd
import csv
from math import cos, pi,sin,sqrt,exp
from dolfin import UserExpression
from dolfin import SpatialCoordinate
from math import cos as mcos, pi as mpi, sin, sqrt, exp 
import ufl                                               



''' Finite-Element-Simulation of a free boundary Stokes flow describing a viscous compressible actin gel
    polymerizing on a curved membrane.

    Copyright: Dennis Woerthmueller, Kristiana Mihali 
    Date last modified: May 12, 2026
'''



# parameters['reorder_dofs_serial'] = True
parameters['allow_extrapolation'] = False
# set some dolfin specific parameters
parameters["form_compiler"]["representation"]="uflacs"
parameters["form_compiler"]["optimize"]=True
parameters["form_compiler"]["cpp_optimize"]=True
# parameters["form_compiler"]["quadrature_degree"]=4

#Create all directories necessary for a single simulation run
# path_to_meshfiles = 'Meshfiles/'
# if not os.path.exists(path_to_meshfiles):
#     os.makedirs(path_to_meshfiles)
import tempfile
path_to_meshfiles = tempfile.mkdtemp() + '/'


path_to_rawData = 'Data/'
if not os.path.exists(path_to_rawData):
    os.makedirs(path_to_rawData)

# read csv i.e. the input parameters for this run
df = pd.read_csv('inputParams.csv')
params = df.iloc[0].to_dict()

gmsh.initialize()

vtkfile_vel = File(path_to_rawData + "vel.pvd")
vtkfile_press = File(path_to_rawData + "press.pvd")
vtkfile_dens = File(path_to_rawData + "dens.pvd")
vtkfile_stress = File(path_to_rawData + "stress.pvd")
vtkfile_trac = File(path_to_rawData + "trac.pvd")

xdmf_flag = True



class SquareCompartment(UserExpression):

    def __init__(self,A, degree=1):
        # print("BIS HIER")
        super().__init__()
        self.A =A
    def eval(self, value, x ):
        eps = 1e-5
        if (x[0] >eps and x[0]< self.A-eps):
             value[0] = 1.0
        else:
            value[0] = 0.0
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        eps = 1e-5
        return on_boundary and (x[0] > eps and x[0] < params['CDLX_ND']-eps) and (near(x[1], 0) or near(x[1], P))

        # def value_shape(self):
        #     return (1,)
def calc_curvature(coordinates):
    x = coordinates[:,0]
    y =coordinates[:,1]
    x_t = np.gradient(coordinates[:, 0])
    y_t = np.gradient(coordinates[:, 1])

    vel = np.array([[x_t[i], y_t[i]] for i in range(x_t.size)])
    speed = np.sqrt(x_t * x_t + y_t * y_t)
    tangent = np.array([1 / speed] * 2).transpose() * vel
    ss_t = np.gradient(speed)
    xx_t = np.gradient(x_t)
    yy_t = np.gradient(y_t)

    C = np.abs(xx_t * y_t - x_t * yy_t) / (x_t * x_t + y_t * y_t) ** 1.5
    return C
def dof_coordinates_sorted(mesh,pbc,boundaries,idx):
    # VDOF = FunctionSpace(mesh, "CG", 1,constrained_domain=pbc)
    VDOF = FunctionSpace(mesh, "CG", 1)
    
    v2d = vertex_to_dof_map(VDOF)
    dofs_boundary = []
    vertices_boundary = []
    for facet in facets(mesh):
        if boundaries[facet.index()] == idx:
            vertices = facet.entities(0)
            for vertex in vertices:
                dofs_boundary.append(v2d[vertex])
                vertices_boundary.append(vertex)
    unique_dofs_boundary = np.array(list(set(dofs_boundary)), dtype=np.int32)
    unique_vertices_boundary= np.array(list(set(vertices_boundary)), dtype=np.int32)
    boundary_coords_boundary = VDOF.tabulate_dof_coordinates()[unique_dofs_boundary]
        # boundary_coords.sort(axis=0)
    # column by which we need to sort the array
    col = 0
    # Sorting by column
    boundary_coords_boundary_sorted = boundary_coords_boundary[np.argsort(boundary_coords_boundary[:,col])]
    sorted_unique_dofs_boundary = unique_dofs_boundary[np.argsort(boundary_coords_boundary[:,col])]
    sorted_unique_vertices_boundary = unique_vertices_boundary[np.argsort(boundary_coords_boundary[:,col])]


    return sorted_unique_dofs_boundary,sorted_unique_vertices_boundary,boundary_coords_boundary_sorted

def calc_normal_cg1_bottom(mesh,ds):
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(-n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L)
    return nh


def calc_tangent_cg1_bottom(mesh,ds):
    n = FacetNormal(mesh)
    t = as_vector([n[1], -n[0]])
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(t, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    th = Function(V)
    solve(A, th.vector(), L)
    return th


def calc_normal_cg2_bottom(mesh,ds):
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(-n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L)
    return nh

def calc_tangent_cg2_bottom(mesh,ds):
    n = FacetNormal(mesh)
    t = as_vector([n[1], -n[0]])
    V = VectorFunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(t, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    th = Function(V)
    solve(A, th.vector(), L)
    return th

def calc_normal_cg2(mesh):
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(-n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L)
    return nh

def calc_tangent_cg2(mesh):
    n = FacetNormal(mesh)
    t = as_vector([n[1], -n[0]])
    V = VectorFunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(t, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    th = Function(V)
    solve(A, nh.vector(), L)
    return th
    

def calc_normal_cg1(mesh):
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(-n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L)
    return nh

def calc_normal_dg0(mesh):
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "DG", 0)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(-n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L)
    return nh

def initial_mesh():

  
    # Constants for domain geommetry
    CDLX = params['CDLX']
    CDLY = params['CDLY']
    delta = params['delta']
    gridsize = params['gridsize']
    q = params['q']
    Amp = params['amp']
    eta = Constant(params['eta'])
    c = params['c'] # only relevant for gauss

    
    
    nPoints = 3000#  Number of discretization points (top-right point of the inlet region)
    P1 = gmsh.model.geo.addPoint(0,Amp,0,gridsize) # cos
    # P1 = gmsh.model.geo.addPoint(0,0,0,gridsize) # gauss
    pointList = [P1]
    for i in np.arange(1,nPoints):
        x = CDLX*i/(nPoints)
        _p = gmsh.model.geo.addPoint(x,Amp*cos(q*x),0,gridsize) #cos 
        # _p = gmsh.model.geo.addPoint(x,-Amp*exp(-pow((x-1),2)/(pow(c,2))),0,gridsize) # gauss
        pointList.append(_p)
    Pend = gmsh.model.geo.addPoint(CDLX,Amp,0,gridsize) # cos
    # Pend = gmsh.model.geo.addPoint(CDLX,0,0,gridsize) # gauss
    pointList.append(Pend)

    
    P1_free = gmsh.model.geo.addPoint(0,Amp+delta,0,gridsize) # cos
    # P1_free = gmsh.model.geo.addPoint(0,delta,0,gridsize) # gauss
    pointList_free = [P1_free]
    for i in np.arange(1,nPoints):
        x = CDLX*i/(nPoints)
        _p = gmsh.model.geo.addPoint(x,Amp*cos(q*x)+delta,0,gridsize) # cos
        # _p = gmsh.model.geo.addPoint(x,-Amp*exp(-pow((x-1),2)/(pow(c,2)))+delta,0,gridsize) # gauss
        pointList_free.append(_p)
    Pend_free = gmsh.model.geo.addPoint(CDLX,Amp+delta,0,gridsize) # cos 
    # Pend_free = gmsh.model.geo.addPoint(CDLX,delta,0,gridsize) #gauss
    pointList_free.append(Pend_free)

    Line1 = gmsh.model.geo.addSpline(pointList)
    Line2 = gmsh.model.geo.addLine(P1, P1_free)
    Line3 = gmsh.model.geo.addSpline(pointList_free)
    Line4 = gmsh.model.geo.addLine(Pend_free, Pend)

    c1 = gmsh.model.geo.addCurveLoop([Line2,Line3,Line4,-Line1])
    planeSurf1 =  gmsh.model.geo.addPlaneSurface([c1])

    physicalSurf1 = gmsh.model.geo.addPhysicalGroup(2, [planeSurf1], tag=1) 
    physicalCurve1 = gmsh.model.geo.addPhysicalGroup(1, [Line1], tag=1) 
    physicalCurve2 = gmsh.model.geo.addPhysicalGroup(1, [Line3], tag=2) 
    physicalCurve3 = gmsh.model.geo.addPhysicalGroup(1, [Line2,Line4], tag=3) 

    gmsh.model.geo.synchronize()
    
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [1])
    gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", 1000)
    
    
    
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", gridsize/10)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", gridsize)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.5)
    
    
    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(7)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)




  
    # gmsh.model.mesh.setTransfiniteCurve(Line2, int(50))
    # gmsh.model.mesh.setTransfiniteCurve(Line4, int(50))
    # gmsh.model.mesh.setTransfiniteCurve(Line1, int(200))
    # gmsh.model.mesh.setTransfiniteCurve(Line3, int(200))
    # gmsh.model.mesh.algorithm = 6
    # gmsh.model.mesh.CharacteristicLengthFromCurvature = 0
    # gmsh.model.mesh.CharacteristicLengthFromPoints = 1

    gmsh.model.mesh.generate(dim=2)
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
    gmsh.write(path_to_meshfiles + "mesh_0.msh")
    gmsh.clear()
    return None

def new_mesh(old_mesh, iteration_step):
    boundaries = MeshFunction("size_t", old_mesh, path_to_meshfiles + "domain_facet_region.xml")
    subdomains = MeshFunction("size_t", old_mesh, path_to_meshfiles + "domain_physical_region.xml")
    V = FunctionSpace(old_mesh, "CG", 1)
    v2d = vertex_to_dof_map(V)

    # save mesh boundary of free surface (top)
    dofs = []
    for facet in facets(old_mesh):
        if boundaries[facet.index()] == 2:
            vertices = facet.entities(0)
            for vertex in vertices:
                dofs.append(v2d[vertex])
    unique_dofs = np.array(list(set(dofs)), dtype=np.int32)
    boundary_coords = V.tabulate_dof_coordinates()[unique_dofs]

    # save mesh boundary of solid surface (bottom)
    dofs_bottom = []
    for facet in facets(old_mesh):
        if boundaries[facet.index()] == 1:
            vertices = facet.entities(0)
            for vertex in vertices:
                dofs_bottom.append(v2d[vertex])
    unique_dofs_bottom = np.array(list(set(dofs_bottom)), dtype=np.int32)
    boundary_coords_bottom = V.tabulate_dof_coordinates()[unique_dofs_bottom]
   

    
    # boundary_coords.sort(axis=0)
    # column by which we need to sort the array
    col = 0
    # Sorting by column
    boundary_coords_sorted = boundary_coords[np.argsort(boundary_coords[:,col])]
    x_vals = boundary_coords_sorted[:,0]
    y_vals = boundary_coords_sorted[:,1]
    boundary_coords_bottom_sorted = boundary_coords_bottom[np.argsort(boundary_coords_bottom[:,col])]
    x_vals_bottom = boundary_coords_bottom_sorted[:,0]
    y_vals_bottom = boundary_coords_bottom_sorted[:,1]

    # Constants for domain geommetry
    CDLX = params['CDLX']
    CDLY = params['CDLY']
    gridsize = params['gridsize']
    q = params['q']
    Amp = params['amp']
    eta = Constant(params['eta'])
    c = params['c'] # only relevant for gauss



    P1 = gmsh.model.geo.addPoint(0,Amp,0,gridsize) # cos
    # P1 = gmsh.model.geo.addPoint(0,0,0,gridsize) # gauss
    pointList = [P1]
    nPoints = 3000#  Number of discretization points (top-right point of the inlet region)
    for i in np.arange(1,nPoints):
        x = CDLX*i/(nPoints)
        _p = gmsh.model.geo.addPoint(x,Amp*cos(q*x),0,gridsize) # cos
        # _p = gmsh.model.geo.addPoint(x,-Amp*exp(-pow((x-1),2)/(pow(c,2))),0,gridsize) # gauss
        pointList.append(_p)

    Pend = gmsh.model.geo.addPoint(CDLX,Amp,0,gridsize) # cos
    # Pend = gmsh.model.geo.addPoint(CDLX,0,0,gridsize) # gauss
    pointList.append(Pend)
    pointList_upperBoundary = []
    for j in range(len(x_vals)):
        _p = gmsh.model.geo.addPoint(x_vals[j],y_vals[j],0,gridsize)
        pointList_upperBoundary.append(_p)

    Line1 = gmsh.model.geo.addSpline(pointList)
    Line2 = gmsh.model.geo.addLine(P1, pointList_upperBoundary[0])
    Line3 = gmsh.model.geo.addSpline(pointList_upperBoundary)
    Line4 = gmsh.model.geo.addLine(pointList_upperBoundary[-1], pointList[-1])

    c1 = gmsh.model.geo.addCurveLoop([Line2,Line3,Line4,-Line1])
    planeSurf1 =  gmsh.model.geo.addPlaneSurface([c1])

    physicalSurf1 = gmsh.model.geo.addPhysicalGroup(2, [planeSurf1], tag=1) 
    physicalCurve1 = gmsh.model.geo.addPhysicalGroup(1, [Line1], tag=1) 
    physicalCurve2 = gmsh.model.geo.addPhysicalGroup(1, [Line3], tag=2) 
    physicalCurve3 = gmsh.model.geo.addPhysicalGroup(1, [Line2,Line4], tag=3) 
    gmsh.model.geo.synchronize()

    # gmsh.model.mesh.setTransfiniteCurve(Line2, int(50))#int(round(2*(y_vals[0]-Amp)/gridsize))
    # gmsh.model.mesh.setTransfiniteCurve(Line4, int(50))
    # gmsh.model.mesh.setTransfiniteCurve(Line2, int(200))
    # gmsh.model.mesh.setTransfiniteCurve(Line4, int(200))
    # gmsh.model.mesh.setTransfiniteCurve(Line1, int(200))            values[0]=0.0 #the other direction is zero

    # gmsh.model.mesh.setTransfiniteCurve(Line3, int(200))
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [1])
    gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", 1000)
    
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", gridsize/10)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", gridsize)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.2)
    
    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(7)
    
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


    gmsh.model.mesh.generate(dim=2)
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
    gmsh.write(path_to_meshfiles + "mesh_%s.msh"%(iteration_step))
    gmsh.clear()

    
    # outputDict = {'x_top':np.array(x_vals),'y_top':np.array(y_vals),'x_bottom':np.array(x_vals_bottom), 'y_bottom':np.array(y_vals_bottom)}
    np.savez(path_to_rawData + 'surfaceProfiles.npz', x_top=np.array(x_vals), y_top=np.array(y_vals),x_bottom=np.array(x_vals_bottom),y_bottom=np.array(y_vals_bottom))
    # df = pd.DataFrame.from_dict(data=outputDict,orient = 'index',dtype='float')
    # df = df.transpose()
    # df = df.astype(float)
    # df.to_csv(path_to_rawData + 'surfaceProfiles.csv',header=True)
    return None



def sim(mesh,iteration_step,dt,remesh=False, save=False):
    
    if remesh is True:
        print("REMESH OF DOMAIN")
        new_mesh(mesh,iteration_step)

        os.system('dolfin-convert %smesh_%i.msh %sdomain.xml'%(path_to_meshfiles,iteration_step,path_to_meshfiles))
        # read mesh and corresoponding bulk and surface subdomains
        mesh = Mesh(path_to_meshfiles + "domain.xml")
    
   
    # define UFL integration measures
    boundaries = MeshFunction("size_t", mesh, path_to_meshfiles + "domain_facet_region.xml")
    subdomains = MeshFunction("size_t", mesh, path_to_meshfiles + "domain_physical_region.xml")
     # Save subdomains
    file_results = XDMFFile(path_to_meshfiles + "subdomains.xdmf")
    file_results.write(subdomains)
    # Save boundaries
    file_results = XDMFFile(path_to_meshfiles + "boundaries.xdmf")
    file_results.write(boundaries)

    dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
    ds_all = Measure('ds', subdomain_data=boundaries)

    # select the no-flux interfaces 
    ds_bottom = ds_all(1)# 
    ds_top = ds_all(2)
    ds_lr = ds_all(3)

    n = FacetNormal(mesh)
    t = as_vector([n[1], -n[0]])
    
    class PeriodicBoundary(SubDomain):

        # Left boundary is "target domain" G
        def inside(self, x, on_boundary):
            return near(x[0], 0) and on_boundary

        # Map right boundary (H) to left boundary (G)
        def map(self, x, y):
            y[0] = x[0] - params['CDLX']
            y[1] = x[1]

    # Create periodic boundary condition
    pbc = PeriodicBoundary()

    #Function spaces for saving
    dFE_DG0 = FiniteElement("DG", mesh.ufl_cell(), 0)  # 0
    dFE_DG1 = FiniteElement("DG", mesh.ufl_cell(), 1)
    dFE_CG1 = FiniteElement("CG", mesh.ufl_cell(), 1)#0
    dFE_CG2 = FiniteElement("CG", mesh.ufl_cell(), 2)
    TensorSpace_DG0 = TensorFunctionSpace(mesh, "DG", 0,constrained_domain=pbc)#0
    TensorSpace_DG1 = TensorFunctionSpace(mesh, "DG", 1,constrained_domain=pbc)  # 0
    TensorSpace_CG1 = TensorFunctionSpace(mesh, "CG", 1,constrained_domain=pbc)  # 0
    TensorSpace_CG2 = TensorFunctionSpace(mesh, "CG", 2,constrained_domain=pbc)  # 0

    # tFE_cg = TensorElement(dFE_cg)
    W_DG0 = TensorSpace_DG0  #FunctionSpace(mesh, TensorSpace)
    W_DG1 = TensorSpace_DG1  # FunctionSpace(mesh, TensorSpace)
    W_CG1 = TensorSpace_CG1  # FunctionSpace(mesh, TensorSpace)
    W_CG2 = TensorSpace_CG2  # FunctionSpace(mesh, TensorSpace)

    K_DG0 = FunctionSpace(mesh, dFE_DG0,constrained_domain=pbc)
    K_DG1 = FunctionSpace(mesh, dFE_DG1,constrained_domain=pbc)
    K_CG1 = FunctionSpace(mesh, dFE_CG1,constrained_domain=pbc)
    K_CG2 = FunctionSpace(mesh, dFE_CG2,constrained_domain=pbc)

    V_CG1 = VectorFunctionSpace(mesh,"CG",1,constrained_domain=pbc)
    V_CG2 = VectorFunctionSpace(mesh,"CG",2,constrained_domain=pbc)

    vel = Function(V_CG2, name='Velocity')
    press = Function(K_CG1, name='Pressure')
    stress = Function(W_CG2, name='Total Stress')
    dens = Function(K_CG1, name='Density')
    trac = Function(V_CG2, name='Traction')
    sig_normal = Function(K_CG2, name='Sig nn')
    sig_tangent = Function(K_CG2, name='Sig nt')
    # stress = Function(W_DG0, name='Total Stress')

    class InflowExpression(UserExpression): #
        def __init__(self, vp, alpha, rho_in, amp, q, **kwargs):
            super().__init__(**kwargs)
            self.vp = vp 
            self.alpha=alpha
            self.rho_in = rho_in
            self.amp=amp
            self.q=q

        def eval(self, values, x):
            v_mod = self.vp - self.alpha*(float(self.q))**2*(self.amp)*cos(float(self.q)*x[0])
            values[0] = self.rho_in * v_mod   # scalar j_in(x)

            
        def value_shape(self):
            return()

    def fixed_boundary_problem(mesh, iteration_step):
        # --- Params (as you had) ---
        vp         = params['vp']
        lmda       = params['lmda']
        kd         = Constant(params['K'])
        amp        = params['amp']
        # n = Constant(params['n'])
        c          = params['c']          # only relevant for gauss (unused here)
        alpha        = Constant(params['alpha'])
        xi   = Constant(params['xi'])
        chi      = Constant(params['chi'])
        q       = Constant(params['q'])
        rho_target = Constant(params['rho_target'])
        rho_in     = Constant(params['rho_in'])
        zeta       = Constant(params['zeta'])
        D          = Constant(params['D'])

        # Optional (set in params if you want): tiny grad-div, friction floor
        # gamma_div_val = float(params.get('gamma_div', 0.0))  # e.g. 0.05–0.2 if ζ=0
        # friction_floor = float(params.get('friction_floor', 1e-8))
        # friction_eff = Constant(max(friction, friction_floor))

        # ---------- Spaces ----------
        V = VectorElement("CG", triangle, 2)
        Q = FiniteElement("CG", triangle, 1)
        W = FunctionSpace(mesh, MixedElement([V, Q]), constrained_domain=pbc)

        U  = Function(W)      # (v, rho)
        Un = Function(W)
        v,  p  = split(U)
        vn, rhon = split(Un)
        w,  ws   = TestFunctions(W)
        I = Identity(mesh.geometry().dim())
        n = FacetNormal(mesh)

        # ---------- Geometry: n,t on the BOTTOM boundary ----------
        # analytical expressions (unit vectors analytically)
        n_exp = Expression(
            ('amp*q*sin(q*x[0])/sqrt(1+pow(amp*q*sin(q*x[0]),2))',
            '1.0/sqrt(1+pow(amp*q*sin(q*x[0]),2))'),
            vp=vp,alpha=alpha, q=q, amp=amp, degree=2
        )
        t_exp = Expression(
            ('1.0/sqrt(1+pow(amp*q*sin(q*x[0]),2))',
            '-amp*q*sin(q*x[0])/sqrt(1+pow(amp*q*sin(q*x[0]),2))'),
            vp=vp, alpha=alpha, q=q, amp=amp, degree=2
        )

        # project to FE space used in assembly & re-orthonormalize AFTER projection
        n_exp_cg2 = project(n_exp, V_CG2)
        t_exp_cg2 = project(t_exp, V_CG2)

        # n_exp_cg2 = project(n_exp_cg2 / sqrt(dot(n_exp_cg2, n_exp_cg2)), V_CG2)
        # t_exp_cg2 = project(t_exp_cg2 - dot(t_exp_cg2, n_exp_cg2) * n_exp_cg2, V_CG2)
        # t_exp_cg2 = project(t_exp_cg2 / sqrt(dot(t_exp_cg2, t_exp_cg2)), V_CG2)

        # Convenience scalar comps (used consistently in assembly)
        vn_sc = dot(v, n_exp_cg2)
        vt_sc = dot(v, t_exp_cg2)
        wn_sc = dot(w, n_exp_cg2)
        wt_sc = dot(w, t_exp_cg2)

        u0 = interpolate(Constant((0.0, 0.0)), V_CG2)
        p0 = interpolate(Constant(0), K_CG1)
        # p0 = interpolate(Constant(0), K_CG1)

        assigner = FunctionAssigner(W,[V_CG2,K_CG1])
        assigner.assign(U, [u0, p0])    

        # ---------- Pressure & stress ----------
        rho = rho_target+p/chi
        # p = alpha * rho

        def stress_function(vv, pp):
            # ζ is your (bulk) dilatational viscosity; ζ=0 allowed
            return -pp*I + zeta*div(vv)*I + sym(grad(vv))

        sigmaU = stress_function(v, p)
        sigmaW = stress_function(w, 0)  # test stress; no pressure part

        # ---------- Nitsche params ----------
        hK     = CellDiameter(mesh)
        gammaN = Constant(80.0)  # your choice; keep 20–200 typical
        v_mod = Expression(
            "vp - alpha*pow(q,2)*amp*cos(q*x[0])",
            vp=float(params["vp"]),
            alpha=float(params["alpha"]),
            q=float(params["q"]),
            amp=float(params["amp"]),
            degree=2,
        )

        # ---------- Momentum weak form (signs kept exactly like your snippet) ----------
        # Natural tractions (split n/t), Robin in t, Nitsche sym+penalty in n
        F_mom = (
            - inner(sigmaU, grad(w)) * dx
            # natural traction components on bottom (use SAME n,t)
            - dot(dot(sigmaU, n_exp_cg2), n_exp_cg2) * wn_sc * ds_bottom
            # - dot(dot(sigmaU, n_exp_cg2), t_exp_cg2) * wt_sc * ds_bottom
            # Robin (tangential) with tiny floor to avoid rigid-mode singularity
            - xi * vt_sc * wt_sc * ds_bottom
            # Nitsche: symmetric consistency + penalty for vn = vp
            # - dot(dot(sigmaW, n_exp_cg2), n_exp_cg2) * (vn_sc - vp) * ds_bottom
            + (gammaN/hK) * (vn_sc - v_mod) * wn_sc * ds_bottom
        )

        # Optional grad-div stabilization if zeta=0 and you see volumetric noise
        # if gamma_div_val > 0.0:
        #     F_mom += Constant(gamma_div_val) * div(v) * div(w) * dx

        # ---------- Mass conservation (conservative) ----------
        p_in = Constant(0.0)
        bc_bottom = DirichletBC(W.sub(1), p_in, boundaries, 1)

        F_p = (
            + kd * p * ws * dx                              
            + chi * kd * rho_target * ws * dx               
            - chi * inner(rho_target * v, grad(ws)) * dx    
            - inner(p * v, grad(ws)) * dx                   
            + D * inner(grad(p), grad(ws)) * dx             
            + (- D * dot(grad(p), n)) * ws * ds_top
            + chi * (dot(rho_target * v, n)) * ws * ds_top 
        )       
        F = F_mom + F_p

        # ---------- Solve ----------
        solve(F == 0, U,[bc_bottom])  
        rho_computed = project(rho, K_CG1)

        # ---------- Diagnostics / post-proc (robust, frame-consistent) ----------
        # Stress & traction in SAME frame used in assembly
        stress_components = project(sigmaU, W_CG1)
        traction_vec = dot(sigmaU, n_exp_cg2)

        # Tangential traction and Robin residual on bottom, projected to DG for clean sampling
        # Vb = FunctionSpace(mesh, "DG", 0)
        # sigma_nt_expr = dot(traction_vec, t_exp_cg2)                   # (σ n)·t
        # sigma_nt_dg   = project(sigma_nt_expr, Vb)
        # vt_dg         = project(vt_sc, Vb)
        # robin_res_dg  = project(sigma_nt_expr + friction*vt_sc, Vb)

        # # (optional) print RMS values
        # Lb = assemble(1.0 * ds_bottom)
        # Sn = np.sqrt(assemble(sigma_nt_expr**2 * ds_bottom)) / np.sqrt(Lb)
        # Rt = np.sqrt(assemble(robin_res_dg**2 * ds_bottom)) / np.sqrt(Lb)
        # print("RMS sigma_nt =", Sn, " | RMS Robin residual =", Rt)

        # ---------- (Your existing exports; kept as-is) ----------
        # Components for CSV (you had CG1/vertex sampling; we keep but also write DG)
        sig_xx = project(stress_components[0,0],  K_CG1)
        sig_xy = project(stress_components[1,0],  K_CG1)
        sig_yy = project(stress_components[1,1],  K_CG1)

        unique_dofs_bottom_sorted, sorted_unique_vertices_bottom, boundary_coords_bottom_sorted = \
            dof_coordinates_sorted(mesh, pbc, boundaries, 1)

        sig_xx_vals = []
        sig_xy_vals = []
        sig_yy_vals = []
        t_x_vals = []
        t_y_vals = []
        n_x_vals = []
        n_y_vals = []
        f_prime_vals = []

        def f_prime(x):
            return - amp * float(q) *np.sin(float(q) * x)

        curvature = calc_curvature(boundary_coords_bottom_sorted)

        for coord in boundary_coords_bottom_sorted:
            x, y = coord
            sig_xx_vals.append(sig_xx((x,y)))
            sig_xy_vals.append(sig_xy((x,y)))
            sig_yy_vals.append(sig_yy((x,y)))
            t_val = t_exp_cg2(x,y)
            n_val = n_exp_cg2(x,y)
            t_x_vals.append(t_val[0]); t_y_vals.append(t_val[1])
            n_x_vals.append(n_val[0]); n_y_vals.append(n_val[1])
            f_prime_vals.append(f_prime(x))

        sig_xx_v = np.flip(np.array(sig_xx_vals))
        sig_xy_v = np.flip(np.array(sig_xy_vals))
        sig_yy_v = np.flip(np.array(sig_yy_vals))
        t_x_v    = np.flip(np.array(t_x_vals))
        t_y_v    = np.flip(np.array(t_y_vals))
        n_x_v    = np.flip(np.array(n_x_vals))
        n_y_v    = np.flip(np.array(n_y_vals))
        f_prime_v= np.flip(np.array(f_prime_vals))

        sig_nn_v = sig_xx_v*n_x_v**2 + 2*sig_xy_v*n_x_v*n_y_v + sig_yy_v*n_y_v**2
        sig_nt_v = sig_xx_v*n_x_v*t_x_v + sig_xy_v*(n_x_v*t_y_v + n_y_v*t_x_v) + sig_yy_v*n_y_v*t_y_v
        f_y_v    = sig_xy_v*n_x_v + sig_yy_v*n_y_v

        # Traction/forces using SAME n,t as assembly
        forceX = (dot(traction_vec, as_vector((1.0,0.0)))) * ds_bottom(degree=7)
        forceY = (dot(traction_vec, as_vector((0.0,1.0)))) * ds_bottom(degree=7)
        forceN = inner(traction_vec, n_exp_cg2) * ds_bottom(degree=7)
        forceT = inner(traction_vec, t_exp_cg2) * ds_bottom(degree=7)
        fX = assemble(forceX); fY = assemble(forceY); fN = assemble(forceN); fT = assemble(forceT)

        # arc-length measure (kept from your code)
        arclength_sin = assemble(inner(dot(Identity(2), n), n) * ds_bottom(degree=7))
        print(arclength_sin)
        print(fX/arclength_sin, fY/arclength_sin, arclength_sin)

        # Export fields (unchanged names)
        vel.assign(project(v, V_CG2))          # nm/s
        press.assign(project(p, K_CG1))        # kPa
        # stress.assign(project(sigmaU, W_CG1))  # kPa
        stress.assign(project(sigmaU, W_CG1))  # kPa
        trac.assign(project(traction_vec, V_CG1))
        dens.assign(project(rho, K_CG1))

        if xdmf_flag is True and save is True:
            xdmf_stokes.write(vel,   iteration_step)
            xdmf_stokes.write(press, iteration_step)
            xdmf_stokes.write(stress,iteration_step)
            xdmf_stokes.write(trac,  iteration_step)
            xdmf_stokes.write(dens,  iteration_step)
            xdmf_stokes.write(n_exp_cg2, iteration_step)  # write the exact frame used

        vtkfile_vel   << (vel,   iteration_step)
        vtkfile_press << (press, iteration_step)
        vtkfile_stress<< (stress,iteration_step)
        vtkfile_dens<< (dens,  iteration_step)
        vtkfile_trac  << (trac,  iteration_step)

        # Save CSV (keep your format) + add robust sigma_nt_dg if you want
        outputDict = {
            'x': boundary_coords_bottom_sorted[:,0], 'y': boundary_coords_bottom_sorted[:,1],
            'sig_xx': sig_xx_v, 'sig_xy': sig_xy_v, 'sig_yy': sig_yy_v,
            'sig_nn': sig_nn_v, 'sig_nt': sig_nt_v, 'f_y': f_y_v,
            'shape_deriv': f_prime_vals, 'C': curvature, 
        }
        pd.DataFrame.from_dict(outputDict).to_csv(path_to_rawData + 'surface_stresses%04d.csv' % (iteration_step), header=True)

        df['fX'] = [fX/arclength_sin]
        df['fY'] = [fY/arclength_sin]
        df.to_csv('./inputParams_add.csv', index=False)

        return U


    
    
    def harmonic_extension(mesh,U_stokes,dt,iteration_step):

        P1 = VectorElement('CG', triangle, 2)
        P2 = FiniteElement('CG',triangle,1)
       
        V_harmonic = FunctionSpace(mesh, P1, constrained_domain=pbc)
        p_harmonic = FunctionSpace(mesh, P2, constrained_domain=pbc)
        # Define functions for variational problem
        du = TrialFunction(V_harmonic)
        u = Function(V_harmonic)

        w= TestFunction(V_harmonic)
        v0 = Constant((0.0,0.0))
        v_stokes = Function(V_CG2, name='Velocity')
        p_stokes = Function(K_CG1, name='Pressure')
        v_stokes.assign(interpolate(U_stokes.sub(0),V_harmonic))
        p_stokes.assign(interpolate(U_stokes.sub(1),p_harmonic))
   
    
        I = Identity(2)
        
        nh = calc_normal_cg2(mesh)
        nh_b = calc_normal_cg2_bottom(mesh,ds_top)
        # bc_sides = DirichletBC(V.sub(0), v0, boundaries, 3) #front
        
        inflow = project(v_stokes,V_harmonic)
        bc_top = DirichletBC(V_harmonic, inflow, boundaries, 2) #front
        bc_bottom = DirichletBC(V_harmonic, v0, boundaries, 1) #front

        
        FWF = inner(grad(u),grad(w))*dx
        solve(FWF==0, u, [bc_top,bc_bottom])
        
        

        
        vel_CG2 = Function(V_harmonic, name='Velocity')
        vel_CG2.assign(u)
     
        if xdmf_flag is True:
            xdmf_harmonic.write(vel_CG2,iteration_step)
            
        ey = Expression(('0.0','1.0'),degree=2)
        uy = dot(u,ey)
        
        
        
        unique_dofs_top_sorted,sorted_unique_vertices_top,boundary_coords_top_sorted = dof_coordinates_sorted(mesh,pbc,boundaries,2)

        vy_vals = []
        p_vals = []
        for coord in boundary_coords_top_sorted:
            x = coord[0]
            y = coord[1]
            vy_vals.append(uy((x,y)))
            p_vals.append(p_stokes((x,y)))


        errMeasure = np.max(np.array(vy_vals))
        n_t = FacetNormal(mesh)
        np.savetxt('errMeasure.txt', [errMeasure])  
        print("ERR MEASURE ", errMeasure)


        
        vel_CG2.assign(project(dt*uy*ey,V_harmonic))
        # un = dot(u,nh_b)
        # vel_CG2.assign(project(dt*un*nh_b,V_harmonic))
        
        
       
        return vel_CG2,errMeasure
    
    
    U_stokes = fixed_boundary_problem(mesh,iteration_step)
    u_harmonic,errMeasure = harmonic_extension(mesh,U_stokes,dt,iteration_step)

    return u_harmonic,mesh,errMeasure


# create an xdmf-file for all the simulation outputs for the Stokes equation
xdmf_stokes = XDMFFile(path_to_rawData + "stokes.xdmf")
xdmf_stokes.parameters["flush_output"] = True
xdmf_stokes.parameters["functions_share_mesh"] = True

# create an xdmf-file for all the simulation outputs for the harmonic problem
xdmf_harmonic = XDMFFile(path_to_rawData + "harmonic.xdmf")
xdmf_harmonic.parameters["flush_output"] = True
xdmf_harmonic.parameters["functions_share_mesh"] = True

initial_mesh()
os.system('dolfin-convert %smesh_0.msh %sdomain.xml'%(path_to_meshfiles,path_to_meshfiles))
# read mesh and corresoponding bulk and surface subdomains
mesh = Mesh(path_to_meshfiles + "domain.xml")

n_step = 0
dt = params['dt'] # for extreme cases this probably should be included in the params csv to be individual for specifications of each run 

current_mesh = mesh
remesh = False

counter = 0
errMeasure = np.inf

for n_step in range(7000):
    if n_step%1 ==0:
        remesh = True

    is_last = (np.absolute(errMeasure) < 1e-5) or (counter > 1000)
    u_harmonic,current_mesh,errMeasure = sim(current_mesh,n_step,dt,remesh, save=is_last)
    print("####################################################### ", errMeasure)
    ALE.move(current_mesh, u_harmonic)
    current_mesh.smooth()
    remesh = False
    if np.absolute(errMeasure) < 1e-5:
        break
        print("SIMULATION FINISHED BY ERR MEASURE")
    elif np.absolute(errMeasure) > 1e-5 and counter > 1000:
        break
        print("SIMULATION FINISHED BY NO OF ITERATION STEPS")
    counter = counter + 1
