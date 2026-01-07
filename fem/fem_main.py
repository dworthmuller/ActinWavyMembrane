from dolfin import *
import numpy as np
import os
import gmsh
import pandas as pd


''' Finite-Element-Simulation of a free boundary Stokes flow describing a viscous actin gel
    polymerizing on a fluctuating (curved) membrane.

    Copyright: Dennis Woerthmueller
    Date last modified: Januray 7, 2026
'''

# parameters['reorder_dofs_serial'] = True  # UNUSED
parameters['allow_extrapolation'] = False
# set some dolfin specific parameters
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True


# -------------------------------------------------------------------
# Create all directories necessary for a single simulation run
# -------------------------------------------------------------------
path_to_meshfiles = 'Meshfiles/'
if not os.path.exists(path_to_meshfiles):
    os.makedirs(path_to_meshfiles)

path_to_rawData = 'Data/'
if not os.path.exists(path_to_rawData):
    os.makedirs(path_to_rawData)

# -------------------------------------------------------------------
# Read csv with input parameters for this run
# -------------------------------------------------------------------
df = pd.read_csv('inputParams.csv')
params = df.iloc[0].to_dict()

# Initialize gmsh once for the whole script
gmsh.initialize()

# VTK output for quick inspection
vtkfile_vel = File(path_to_rawData + "vel.pvd")
vtkfile_press = File(path_to_rawData + "press.pvd")
vtkfile_stress = File(path_to_rawData + "stress.pvd")
vtkfile_trac = File(path_to_rawData + "trac.pvd")

xdmf_flag = True



# -------------------------------------------------------------------
# Geometry / curvature helpers
# -------------------------------------------------------------------

def calc_curvature(coordinates):
    """Compute curvature of a parametric curve given as (x,y) point array."""
    x = coordinates[:, 0]
    y = coordinates[:, 1]
    x_t = np.gradient(x)
    y_t = np.gradient(y)

    vel = np.array([[x_t[i], y_t[i]] for i in range(x_t.size)])
    speed = np.sqrt(x_t * x_t + y_t * y_t)
    tangent = np.array([1 / speed] * 2).transpose() * vel  # UNUSED variable 'tangent'

    ss_t = np.gradient(speed)
    xx_t = np.gradient(x_t)
    yy_t = np.gradient(y_t)

    C = np.abs(xx_t * y_t - x_t * yy_t) / (x_t * x_t + y_t * y_t) ** 1.5
    return C


def dof_coordinates_sorted(mesh, pbc, boundaries, idx):
    """
    Collect and sort boundary DOFs and vertices for a given boundary index 'idx'
    along the x-coordinate.
    """
    # VDOF = FunctionSpace(mesh, "CG", 1, constrained_domain=pbc)
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
    unique_vertices_boundary = np.array(list(set(vertices_boundary)), dtype=np.int32)
    boundary_coords_boundary = VDOF.tabulate_dof_coordinates()[unique_dofs_boundary]

    # Sort by x-coordinate
    col = 0
    boundary_coords_boundary_sorted = boundary_coords_boundary[np.argsort(boundary_coords_boundary[:, col])]
    sorted_unique_dofs_boundary = unique_dofs_boundary[np.argsort(boundary_coords_boundary[:, col])]
    sorted_unique_vertices_boundary = unique_vertices_boundary[np.argsort(boundary_coords_boundary[:, col])]

    return sorted_unique_dofs_boundary, sorted_unique_vertices_boundary, boundary_coords_boundary_sorted


# --- normal/tangent reconstruction helpers on boundaries ---
def calc_normal_cg2_bottom(mesh, ds):
    """Compute CG2 approximation of outward normal on given boundary measure ds."""
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u, v) * ds
    l = inner(-n, v) * ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L)
    return nh


def calc_tangent_cg2_bottom(mesh, ds):
    """Compute CG2 approximation of tangent on given boundary measure ds."""
    n = FacetNormal(mesh)
    t = as_vector([n[1], -n[0]])
    V = VectorFunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u, v) * ds
    l = inner(t, v) * ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    th = Function(V)
    solve(A, th.vector(), L)
    return th


def calc_normal_cg2(mesh):
    """Compute CG2 approximation of normal on all boundaries (via ds)."""
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u, v) * ds
    l = inner(-n, v) * ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L)
    return nh



def calc_normal_cg1(mesh):
    """Compute CG1 approximation of normal on all boundaries (used for output)."""
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u, v) * ds
    l = inner(-n, v) * ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)
    solve(A, nh.vector(), L)
    return nh




# -------------------------------------------------------------------
# Initial mesh construction (cosine or Gaussian top surface)
# -------------------------------------------------------------------

def initial_mesh():
    """
    Build the initial 2D domain mesh with a prescribed bottom and free (top) surface
    using gmsh and save it as 'mesh_0.msh'.
    """

    # Constants for domain geometry
    CDLX = params['CDLX']
    delta = params['delta']
    gridsize = params['gridsize']
    q = params['q']
    Amp = params['amp']
    c = params['c']  # only relevant for gaussian shape (currently commented out)

    nPoints = 3000  # Number of discretization points along the interface

    # Bottom (solid) surface: cosine
    P1 = gmsh.model.geo.addPoint(0, Amp, 0, gridsize)  # cos
    # P1 = gmsh.model.geo.addPoint(0, 0, 0, gridsize)  # gauss
    pointList = [P1]
    for i in np.arange(1, nPoints):
        x = CDLX * i / (nPoints)
        _p = gmsh.model.geo.addPoint(x, Amp * cos(q * x), 0, gridsize)  # cos
        # _p = gmsh.model.geo.addPoint(x, -Amp*exp(-pow((x-1),2)/(pow(c,2))), 0, gridsize)  # gauss
        pointList.append(_p)
    Pend = gmsh.model.geo.addPoint(CDLX, Amp, 0, gridsize)  # cos
    # Pend = gmsh.model.geo.addPoint(CDLX, 0, 0, gridsize)  # gauss
    pointList.append(Pend)

    # Top (free) surface: offset by delta
    P1_free = gmsh.model.geo.addPoint(0, Amp + delta, 0, gridsize)  # cos
    # P1_free = gmsh.model.geo.addPoint(0, delta, 0, gridsize)  # gauss
    pointList_free = [P1_free]
    for i in np.arange(1, nPoints):
        x = CDLX * i / (nPoints)
        _p = gmsh.model.geo.addPoint(x, Amp * cos(q * x) + delta, 0, gridsize)  # cos
        # _p = gmsh.model.geo.addPoint(x, -Amp*exp(-pow((x-1),2)/(pow(c,2)))+delta, 0, gridsize)  # gauss
        pointList_free.append(_p)
    Pend_free = gmsh.model.geo.addPoint(CDLX, Amp + delta, 0, gridsize)  # cos
    # Pend_free = gmsh.model.geo.addPoint(CDLX, delta, 0, gridsize)  # gauss
    pointList_free.append(Pend_free)

    # Build 2D surface
    Line1 = gmsh.model.geo.addSpline(pointList)
    Line2 = gmsh.model.geo.addLine(P1, P1_free)
    Line3 = gmsh.model.geo.addSpline(pointList_free)
    Line4 = gmsh.model.geo.addLine(Pend_free, Pend)

    c1 = gmsh.model.geo.addCurveLoop([Line2, Line3, Line4, -Line1])
    planeSurf1 = gmsh.model.geo.addPlaneSurface([c1])

    # Tag physical regions: 2D domain and 1D boundaries
    physicalSurf1 = gmsh.model.geo.addPhysicalGroup(2, [planeSurf1], tag=1)
    physicalCurve1 = gmsh.model.geo.addPhysicalGroup(1, [Line1], tag=1)  # bottom
    physicalCurve2 = gmsh.model.geo.addPhysicalGroup(1, [Line3], tag=2)  # top (free)
    physicalCurve3 = gmsh.model.geo.addPhysicalGroup(1, [Line2, Line4], tag=3)  # left/right

    gmsh.model.geo.synchronize()

    # Mesh refinement near curve #1 (bottom)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [1])
    gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", 1000)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", gridsize / 50)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", gridsize)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.5)

    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(7)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(dim=2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(path_to_meshfiles + "mesh_0.msh")
    gmsh.clear()
    return None


# -------------------------------------------------------------------
# Remeshing based on updated free surface
# -------------------------------------------------------------------

def new_mesh(old_mesh, iteration_step):
    """
    Build a new gmsh mesh whose top boundary follows the current free surface
    extracted from 'old_mesh', then write 'mesh_<iteration_step>.msh'.
    """
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

    # Sort by x-coordinate
    col = 0
    boundary_coords_sorted = boundary_coords[np.argsort(boundary_coords[:, col])]
    x_vals = boundary_coords_sorted[:, 0]
    y_vals = boundary_coords_sorted[:, 1]
    boundary_coords_bottom_sorted = boundary_coords_bottom[np.argsort(boundary_coords_bottom[:, col])]
    x_vals_bottom = boundary_coords_bottom_sorted[:, 0]
    y_vals_bottom = boundary_coords_bottom_sorted[:, 1]

    # Constants for domain geometry (cosine reference shape)
    CDLX = params['CDLX']
    CDLY = params['CDLY']  # UNUSED here
    gridsize = params['gridsize']
    q = params['q']
    Amp = params['amp']
    c = params['c']  # only relevant for gauss

    P1 = gmsh.model.geo.addPoint(0, Amp, 0, gridsize)  # bottom ref: cos
    # P1 = gmsh.model.geo.addPoint(0, 0, 0, gridsize)  # gauss
    pointList = [P1]
    nPoints = 3000  # Number of discretization points
    for i in np.arange(1, nPoints):
        x = CDLX * i / (nPoints)
        _p = gmsh.model.geo.addPoint(x, Amp * cos(q * x), 0, gridsize)  # cos
        # _p = gmsh.model.geo.addPoint(x, -Amp*exp(-pow((x-1),2)/(pow(c,2))), 0, gridsize)  # gauss
        pointList.append(_p)

    Pend = gmsh.model.geo.addPoint(CDLX, Amp, 0, gridsize)  # cos
    # Pend = gmsh.model.geo.addPoint(CDLX, 0, 0, gridsize)  # gauss
    pointList.append(Pend)

    # New upper boundary following numerical free surface
    pointList_upperBoundary = []
    for j in range(len(x_vals)):
        _p = gmsh.model.geo.addPoint(x_vals[j], y_vals[j], 0, gridsize)
        pointList_upperBoundary.append(_p)

    Line1 = gmsh.model.geo.addSpline(pointList)                 # reference bottom
    Line2 = gmsh.model.geo.addLine(P1, pointList_upperBoundary[0])
    Line3 = gmsh.model.geo.addSpline(pointList_upperBoundary)   # updated top
    Line4 = gmsh.model.geo.addLine(pointList_upperBoundary[-1], pointList[-1])

    c1 = gmsh.model.geo.addCurveLoop([Line2, Line3, Line4, -Line1])
    planeSurf1 = gmsh.model.geo.addPlaneSurface([c1])

    physicalSurf1 = gmsh.model.geo.addPhysicalGroup(2, [planeSurf1], tag=1)
    physicalCurve1 = gmsh.model.geo.addPhysicalGroup(1, [Line1], tag=1)
    physicalCurve2 = gmsh.model.geo.addPhysicalGroup(1, [Line3], tag=2)
    physicalCurve3 = gmsh.model.geo.addPhysicalGroup(1, [Line2, Line4], tag=3)
    gmsh.model.geo.synchronize()

    # Mesh size field (similar idea as initial mesh)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [1])
    gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", 1000)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", gridsize / 50)
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
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(path_to_meshfiles + "mesh_%s.msh" % (iteration_step))
    gmsh.clear()

    # Save current top/bottom profiles for post-processing
    np.savez(path_to_rawData + 'surfaceProfiles.npz',
             x_top=np.array(x_vals), y_top=np.array(y_vals),
             x_bottom=np.array(x_vals_bottom), y_bottom=np.array(y_vals_bottom))

    return None


# -------------------------------------------------------------------
# Single time step: solve Stokes + harmonic extension for ALE update
# -------------------------------------------------------------------

def sim(mesh, iteration_step, dt, remesh=False):

    if remesh is True:
        print("REMESH OF DOMAIN")
        new_mesh(mesh, iteration_step)
        os.system('dolfin-convert %smesh_%i.msh %sdomain.xml' % (path_to_meshfiles, iteration_step, path_to_meshfiles))
        # read mesh and corresponding bulk and surface subdomains
        mesh = Mesh(path_to_meshfiles + "domain.xml")

    # define UFL integration measures and read boundary markers
    boundaries = MeshFunction("size_t", mesh, path_to_meshfiles + "domain_facet_region.xml")
    subdomains = MeshFunction("size_t", mesh, path_to_meshfiles + "domain_physical_region.xml")

    # Save subdomains and boundaries for reference (XDMF)
    file_results = XDMFFile(path_to_meshfiles + "subdomains.xdmf")
    file_results.write(subdomains)
    file_results = XDMFFile(path_to_meshfiles + "boundaries.xdmf")
    file_results.write(boundaries)

    dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
    ds_all = Measure('ds', subdomain_data=boundaries)

    # select the no-flux interfaces
    ds_bottom = ds_all(1)
    ds_top = ds_all(2)

    n = FacetNormal(mesh)


    class PeriodicBoundary(SubDomain):
        """Periodic boundary condition in x-direction."""

        # Left boundary is "target domain" G
        def inside(self, x, on_boundary):
            return near(x[0], 0) and on_boundary

        # Map right boundary (H) to left boundary (G)
        def map(self, x, y):
            y[0] = x[0] - params['CDLX']
            y[1] = x[1]

    # Create periodic boundary condition
    pbc = PeriodicBoundary()

    # Function spaces for saving output
    dFE_CG1 = FiniteElement("CG", mesh.ufl_cell(), 1)
    dFE_CG2 = FiniteElement("CG", mesh.ufl_cell(), 2)

    TensorSpace_CG1 = TensorFunctionSpace(mesh, "CG", 1, constrained_domain=pbc)

    W_CG1 = TensorSpace_CG1
    K_CG1 = FunctionSpace(mesh, dFE_CG1, constrained_domain=pbc)
    K_CG2 = FunctionSpace(mesh, dFE_CG2, constrained_domain=pbc)
    V_CG1 = VectorFunctionSpace(mesh, "CG", 1, constrained_domain=pbc)
    V_CG2 = VectorFunctionSpace(mesh, "CG", 2, constrained_domain=pbc)

    # Output functions
    vel = Function(V_CG2, name='Velocity')
    press = Function(K_CG1, name='Pressure')
    stress = Function(W_CG1, name='Total Stress')
    trac = Function(V_CG2, name='Traction')



    # -----------------------------------------------------------------
    # Stokes problem on a fixed domain + post-processing
    # -----------------------------------------------------------------
    def fixed_boundary_problem(mesh, iteration_step):

        # Mixed finite elements: velocity (P2) + pressure (P1)
        P1 = VectorElement('CG', triangle, 2)
        P2 = FiniteElement('CG', triangle, 1)
        element = MixedElement([P1, P2])

        V = FunctionSpace(mesh, element, constrained_domain=pbc)

        # Define functions for variational problem
        dU = TrialFunction(V)
        U_tot = Function(V)
        w, ws = TestFunctions(V)

        # Initial condition for velocity (used only in BCs)
        v0 = Constant((0.0, 0.0))
        p0 = Constant(0)  # UNUSED

        # Split system functions to access components
        v, p = split(U_tot)
        I = Identity(2)

        # Physical parameters
        vp = params['vp']
        q = params['q']
        alpha = params['alpha']      # scalar
        K = Constant(params['K'])    # incompressibility penalty / source term
        amp = params['amp']
        c = params['c']  # only relevant for gauss (commented)

        # Geometric normals and tangents on boundaries
        nh_b = calc_normal_cg2_bottom(mesh, ds_bottom)
        th_b = calc_tangent_cg2_bottom(mesh, ds_bottom)
        nh = calc_normal_cg1(mesh)


        alpha = Constant(params['alpha'])  # overrides scalar alpha, used in Expressions below

        # --- Analytical normal and tangent on cosine shape (bottom) ---
        # f'(x) = -amp*q*sin(q x)
        def f_prime(x):
            return - amp * q * np.sin(q * x)

        # Unit normal (pointing upward) for cos interface
        n_exp = Expression(
            ('amp*q*sin(q*x[0]) / sqrt(1+pow(amp*q*sin(q*x[0]),2))',
             '1 / sqrt(1+pow(amp*q*sin(q*x[0]),2))'),
            vp=vp, amp=amp, alpha=alpha, q=q, degree=2
        )
        n_exp_cg2 = project(n_exp, V_CG2)

        # Unit tangent for cos interface
        t_exp = Expression(
            ('1 / sqrt(1+pow(amp*q*sin(q*x[0]),2))',
             '-amp*q*sin(q*x[0]) / sqrt(1+pow(amp*q*sin(q*x[0]),2))'),
            p=vp, amp=amp, alpha=alpha, q=q, degree=2
        )
        t_exp_cg2 = project(t_exp, V_CG2)

        # Normal inflow velocity: (vp - alpha*q^2*amp*cos(qx)) * n
        inflow = project(
            Expression('vp - alpha*pow(q,2)*amp*cos(q*x[0])',
                       vp=vp, amp=amp, alpha=alpha, q=q, degree=2) * n_exp,
            V_CG2
        )

        # For gaussian shapes, alternative expressions are commented out
        # ...

        # Dirichlet BC: bottom velocity prescribed as 'inflow'
        bc_bottom = DirichletBC(V.sub(0), inflow, boundaries, 1)

        # Stokes weak form with simple pressure constraint (div v = -K)
        FWF = -inner(sym(grad(v)), grad(w)) * dx + p * div(w) * dx + ws * div(v) * dx - (-K) * ws * dx

        # Note: only bc_bottom is currently applied (bc_sides omitted)
        solve(FWF == 0, U_tot, [bc_bottom])

        # Total Cauchy stress (viscous + isotropic pressure)
        total_stress = -p * I + sym(grad(v))

        # Project stress to CG1 tensor space and extract components
        stress_components = project(total_stress, W_CG1)
        sig_xx = project(stress_components[0, 0], K_CG1)
        sig_yy = project(stress_components[1, 1], K_CG1)
        sig_xy = project(stress_components[1, 0], K_CG1)

        # Bottom boundary coordinates, sorted along x
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

        # Curvature of bottom boundary (numeric)
        curvature = calc_curvature(boundary_coords_bottom_sorted)

        # Sample stresses and analytic normal/tangent on bottom boundary
        for coord in boundary_coords_bottom_sorted:
            x = coord[0]
            y = coord[1]
            sig_xx_vals.append(sig_xx((x, y)))
            sig_xy_vals.append(sig_xy((x, y)))
            sig_yy_vals.append(sig_yy((x, y)))
            t_val = t_exp_cg2(x, y)   # analytic tangent
            n_val = n_exp_cg2(x, y)   # analytic normal
            t_x_vals.append(t_val[0])
            t_y_vals.append(t_val[1])
            n_x_vals.append(n_val[0])
            n_y_vals.append(n_val[1])
            f_prime_vals.append(f_prime(x))

        # Flip arrays (to use consistent left-right orientation)
        sig_xx_v = np.flip(np.array(sig_xx_vals))
        sig_xy_v = np.flip(np.array(sig_xy_vals))
        sig_yy_v = np.flip(np.array(sig_yy_vals))
        t_x_v = np.flip(np.array(t_x_vals))
        t_y_v = np.flip(np.array(t_y_vals))
        n_x_v = np.flip(np.array(n_x_vals))
        n_y_v = np.flip(np.array(n_y_vals))
        f_prime_v = np.flip(np.array(f_prime_vals))

        # Normal and tangential stress components (scalar)
        sig_nn_v = sig_xx_v * n_x_v**2 + 2 * sig_xy_v * n_x_v * n_y_v + sig_yy_v * n_y_v**2
        sig_nt_v = sig_xx_v * n_x_v * t_x_v + sig_xy_v * (n_x_v * t_y_v + n_y_v * t_x_v) + sig_yy_v * n_y_v * t_y_v

        # Force in y from stress tensor and normal
        f_y_v = sig_xy_v * n_x_v + sig_yy_v * n_y_v

        # Traction vector on bottom boundary
        traction = dot(total_stress, nh_b)
        forceX = traction[0] * ds_bottom(degree=7)
        forceY = traction[1] * ds_bottom(degree=7)
        forceN = inner(traction, nh_b) * ds_bottom(degree=7)
        forceT = inner(traction, th_b) * ds_bottom(degree=7)
        fX = assemble(forceX)
        fY = assemble(forceY)
        fN = assemble(forceN)
        fT = assemble(forceT)


        # Effective arclength measure (integral of |n|^2 = 1 over ds gives length)
        arclength_sin = assemble(inner(dot(Identity(2), n), n) * ds_bottom(degree=7))
        print(arclength_sin)
        print(fX / arclength_sin, fY / arclength_sin, arclength_sin)

        # Output fields
        vel.assign(project(v, V_CG2))      
        press.assign(project(p, K_CG1))    
        stress.assign(project(total_stress, W_CG1))  
        trac.assign(project(traction, V_CG1))

        if xdmf_flag is True:
            xdmf_stokes.write(vel, iteration_step)
            xdmf_stokes.write(press, iteration_step)
            xdmf_stokes.write(stress, iteration_step)
            xdmf_stokes.write(trac, iteration_step)
            xdmf_stokes.write(nh, iteration_step)
            xdmf_stokes.write(nh_b, iteration_step)

        vtkfile_vel << (vel, iteration_step)
        vtkfile_press << (press, iteration_step)
        vtkfile_stress << (stress, iteration_step)
        vtkfile_trac << (trac, iteration_step)

        # Save boundary-resolved stress data along bottom interface
        outputDict = {
            'x': boundary_coords_bottom_sorted[:, 0],
            'y': boundary_coords_bottom_sorted[:, 1],
            'sig_xx': sig_xx_v,
            'sig_xy': sig_xy_v,
            'sig_yy': sig_yy_v,
            'sig_nn': sig_nn_v,
            'sig_nt': sig_nt_v,
            'f_y': f_y_v,
            'shape_deriv': f_prime_vals,
            'C': curvature
        }
        pd.DataFrame.from_dict(data=outputDict).to_csv(
            path_to_rawData + 'surface_stresses%04d.csv' % (iteration_step),
            header=True
        )

        # Add average forces back to the parameter table (overwrites df)
        df['fX'] = [fX / arclength_sin]
        df['fY'] = [fY / arclength_sin]
        df.to_csv('./inputParams_add.csv', index=False)

        return U_tot

    # -----------------------------------------------------------------
    # Harmonic extension of boundary displacement (ALE mesh update)
    # -----------------------------------------------------------------
    def harmonic_extension(mesh, U_stokes, dt, iteration_step):
        """
        Compute harmonic extension of the Stokes velocity field from the top boundary
        into the bulk, with zero displacement at the bottom. Returns dt * normal
        displacement as ALE update field and an error measure (max vy at top).
        """

        P1 = VectorElement('CG', triangle, 2)
        V_harmonic = FunctionSpace(mesh, P1, constrained_domain=pbc)

        u = Function(V_harmonic)
        w = TestFunction(V_harmonic)

        v0 = Constant((0.0, 0.0))
        v_stokes = Function(V_CG2, name='Velocity')
        v_stokes.assign(interpolate(U_stokes.sub(0), V_harmonic))

        I = Identity(2)  # UNUSED
        nh = calc_normal_cg2(mesh)          # UNUSED
        nh_t = calc_normal_cg2_bottom(mesh, ds_top)

        # Harmonic extension: u solves Laplace(u) = 0 with Dirichlet BCs
        inflow = project(v_stokes, V_harmonic)
        bc_top = DirichletBC(V_harmonic, inflow, boundaries, 2)
        bc_bottom = DirichletBC(V_harmonic, v0, boundaries, 1)

        FWF = inner(grad(u), grad(w)) * dx
        solve(FWF == 0, u, [bc_top, bc_bottom])

        vel_CG2 = Function(V_harmonic, name='Velocity')
        vel_CG2.assign(u)

        # Error measure: maximum vertical velocity at the top boundary
        ey = Expression(('0.0', '1.0'), degree=2)
        uy = dot(u, ey)
        un = dot(u, nh_t)

        unique_dofs_top_sorted, sorted_unique_vertices_top, boundary_coords_top_sorted = \
            dof_coordinates_sorted(mesh, pbc, boundaries, 2)

        vy_vals = []
        for coord in boundary_coords_top_sorted:
            x = coord[0]
            y = coord[1]
            vy_vals.append(uy((x, y)))

        errMeasure = np.max(vy_vals)
        np.savetxt('errMeasure.txt', [errMeasure])
        print("ERR MEASURE ", errMeasure)

        # ALE displacement: only normal component times dt
        vel_CG2.assign(project(dt * un * ey, V_harmonic))
        # Alternative: dt * un * nh_b

        if xdmf_flag is True:
            xdmf_harmonic.write(vel_CG2, iteration_step)

        return vel_CG2, errMeasure

    # ---- run Stokes + harmonic extension ----
    U_stokes = fixed_boundary_problem(mesh, iteration_step)
    u_harmonic, errMeasure = harmonic_extension(mesh, U_stokes, dt, iteration_step)

    return u_harmonic, mesh, errMeasure


# -------------------------------------------------------------------
# XDMF files for time series output
# -------------------------------------------------------------------
xdmf_stokes = XDMFFile(path_to_rawData + "stokes.xdmf")
xdmf_stokes.parameters["flush_output"] = True
xdmf_stokes.parameters["functions_share_mesh"] = True

xdmf_harmonic = XDMFFile(path_to_rawData + "harmonic.xdmf")
xdmf_harmonic.parameters["flush_output"] = True
xdmf_harmonic.parameters["functions_share_mesh"] = True

# -------------------------------------------------------------------
# Main time-stepping loop: remesh + ALE update until convergence
# -------------------------------------------------------------------
initial_mesh()
os.system('dolfin-convert %smesh_0.msh %sdomain.xml' % (path_to_meshfiles, path_to_meshfiles))
mesh = Mesh(path_to_meshfiles + "domain.xml")

n = 0
dt = params['dt']  

current_mesh = mesh
remesh = False

counter = 0
for n in range(7000):
    if n % 1 == 0:   # always true; effectively remesh every iteration
        remesh = True

    u_harmonic, current_mesh, errMeasure = sim(current_mesh, n, dt, remesh)
    print("####################################################### ", errMeasure)

    # ALE move + mesh smoothing
    ALE.move(current_mesh, u_harmonic)
    current_mesh.smooth()
    remesh = False

    # Stopping conditions
    if np.absolute(errMeasure) < 1e-6:
        break
    elif np.absolute(errMeasure) > 1e-6 and counter > 1000:
        break
    counter = counter + 1
