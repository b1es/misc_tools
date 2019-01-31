import vtk
import numpy as np

from vtk.util import numpy_support

from vtk_utils import read_vtk




def project_wss_derivative(du = 0.001, velocity_file=None, stlfile='c0006.stl',  output_fn='surface.vtp',
                           velocity_name = 'v [m/s]',
                           write=False, mu=1.0, move_mesh=False, sign=1):
    '''
    Equidistant points from stl in normal direction
    Takes stl and vtu/vti and generates vtp  
    '''
    stl = read_vtk(stlfile)
    vertices = numpy_support.vtk_to_numpy(stl.GetPoints().GetData())
    indices = numpy_support.vtk_to_numpy(stl.GetPolys().GetData()).reshape(-1, 4)[:, 1:4]
    merged = vtk.vtkPolyData()
    merged.DeepCopy(stl)
    
    vel_data = read_vtk(velocity_file)
    #vfm = getFileReaderOutput(fvmfile)
    

    # Compute normals to vertices
    normalGenerator = vtk.vtkPolyDataNormals()
    normalGenerator.SetInputData(merged)
    normalGenerator.ComputePointNormalsOn()
    normalGenerator.ComputeCellNormalsOff()
    normalGenerator.SetSplitting(0)
    normalGenerator.SetConsistency(0)
    normalGenerator.Update()

    merged = normalGenerator.GetOutput()
    normals = numpy_support.vtk_to_numpy(merged.GetPointData().GetNormals())

    
    # WSS - LBM
    points = vtk.vtkPoints()
    pointsPolyData = vtk.vtkPolyData()

    for normal, pos in zip(normals, vertices):
        points.InsertNextPoint(pos + normal * (-du*sign))

    pointsPolyData.SetPoints(points)
    probe_filter = vtk.vtkProbeFilter()
    probe_filter.SetInputData(pointsPolyData)
    probe_filter.SetSourceData(vel_data)
    probe_filter.GetOutputPort()
    probe_filter.Update()

    probed_data = probe_filter.GetOutput().GetPointData()


    if isinstance( velocity_name, str):
        v_vec = numpy_support.vtk_to_numpy(probed_data.GetArray(velocity_name))
    elif isinstance( velocity_name, list):
        vx = numpy_support.vtk_to_numpy(probed_data.GetArray(velocity_name[0]))
        vy = numpy_support.vtk_to_numpy(probed_data.GetArray(velocity_name[1]))
        vz = numpy_support.vtk_to_numpy(probed_data.GetArray(velocity_name[2]))
        v_vec = np.stack((vx, vy, vz), 1).reshape((-1, 3))
    else:
        raise NotImplementedError
        

    stress_norm = vtk.vtkFloatArray()
    stress_norm.SetNumberOfComponents(1)
    stress_norm.SetName("wss_n")

    
    for v in v_vec:
        s = mu * v/du

        if np.max(s) > 1e35 or np.min(s) < -1e35:
            s = np.array([np.nan])

        stress_norm.InsertNextTypedTuple( [np.linalg.norm(s)])

    merged.GetPointData().AddArray(stress_norm)



    velocity = vtk.vtkFloatArray()
    velocity.SetNumberOfComponents(3)
    velocity.SetName("velocity_at_epsilon")

    
    for v in v_vec:

        if np.max(v) > 1e33 or np.min(v) < -1e33:
            s = np.array( 3*[0.001])

        velocity.InsertNextTypedTuple( v )

    merged.GetPointData().AddArray(velocity)

    if move_mesh:
        merged.SetPoints(points)
        print('moving point by dx inside')
    
    merged.Modified()
    if write:
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(output_fn)
        writer.SetInputData(merged)
        writer.Write()
        print(output_fn," written")

    return merged




