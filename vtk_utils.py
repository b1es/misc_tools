import vtk 
from vtk.util import numpy_support
import numpy as np 

print("VTK version:",vtk.VTK_MAJOR_VERSION)

def vtk2npy(vtkfile,verbose=False,rho_name='rho',v_name='v'):
    """
    Read vti files from sailfish and return 
    gets only rho and v
    """
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(vtkfile)
    reader.Update()
    imageData = reader.GetOutput()
    info = reader.GetInformation()

    field_names = [reader.GetPointArrayName(i) for  i in range(reader.GetNumberOfPointArrays())]
    if verbose:
        print( "reading:",reader.GetFileName() ) 
        print( "fields:",field_names)

    assert( rho_name in field_names and v_name in field_names )

    dims = imageData.GetDimensions()

    data_rho = numpy_support.vtk_to_numpy(imageData.GetPointData().GetArray(rho_name))
    data_v = numpy_support.vtk_to_numpy(imageData.GetPointData().GetArray(v_name))

    data_rho = data_rho.reshape(dims[2], dims[1], dims[0])
    data_v = data_v.reshape(dims[2], dims[1], dims[0],3)
    data_v = np.rollaxis(data_v,3,0)
    # alt.  data = data.transpose(2,1,0)
    
    return (data_rho,data_v)


def probe_vtu(vtu_file='output.vtu',point_data=[[0.050640, 0.027959, 0.05213]],fname=None,shape2d=None,verbose=False):
    '''
    get values of interpolated from vtu mesh on the set of points (N,3)
    '''
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file)
    reader.Update()

    points = vtk.vtkPoints()

    for point in point_data:
        points.InsertNextPoint(*point)

    polydata  = vtk.vtkPolyData()
    polydata.SetPoints(points)

    probe = vtk.vtkProbeFilter()
    probe.SetSourceConnection(reader.GetOutputPort())
    probe.SetInputData(polydata)
    probe.Update()

    out = probe.GetOutput()
    # out.GetBounds()
    ## to get points: numpy_support.vtk_to_numpy(out.GetPoints().GetData()).shape
    pd = out.GetAttributesAsFieldData(0)
    if fname==None:
        # all fields
        output = dict()
        for i in range(pd.GetNumberOfArrays()):
            v_interp_on_grid = numpy_support.vtk_to_numpy(pd.GetArray(i))   
            if shape2d:
                v_interp_on_grid = v_interp_on_grid.reshape(shape2d)
            if verbose:
                print("appending in output:",pd.GetArrayName(i) )
            output[pd.GetArrayName(i)] = v_interp_on_grid
        assert(len(output)>0)    
        return output
    else:    
        field_numbers  = [i for i in range(pd.GetNumberOfArrays()) if pd.GetArrayName(i)==fname]
        assert(len(field_numbers)==1)
        if verbose:
            print("output:",pd.GetArrayName(field_numbers[0]) )
        v_interp_on_grid = numpy_support.vtk_to_numpy(pd.GetArray(field_numbers[0]))   
        if shape2d:
            return v_interp_on_grid.reshape(shape2d)
        return v_interp_on_grid





def probe_vti(vti_file='output.vti',point_data=[[0.050640, 0.027959, 0.05213]],fname=None,shape2d=None,verbose=False):
    '''
    get values of interpolated from vti file mesh on the set of points (N,3)
    '''
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(vti_file)
    reader.Update()

    points = vtk.vtkPoints()

    for point in point_data:
        points.InsertNextPoint(*point)

    polydata  = vtk.vtkPolyData()
    polydata.SetPoints(points)

    probe = vtk.vtkProbeFilter()
    probe.SetSourceConnection(reader.GetOutputPort())
    probe.SetInputData(polydata)
    probe.Update()

    out = probe.GetOutput()
    # out.GetBounds()
    ## to get points: numpy_support.vtk_to_numpy(out.GetPoints().GetData()).shape
    pd = out.GetAttributesAsFieldData(0)
    if fname==None:
        # all fields
        output = dict()
        for i in range(pd.GetNumberOfArrays()):
            v_interp_on_grid = numpy_support.vtk_to_numpy(pd.GetArray(i))   
            if shape2d:
                v_interp_on_grid = v_interp_on_grid.reshape(shape2d)
            if verbose:
                print("appending in output:",pd.GetArrayName(i) )
            output[pd.GetArrayName(i)] = v_interp_on_grid
        assert(len(output)>0)    
        return output
    else:    
        field_numbers  = [i for i in range(pd.GetNumberOfArrays()) if pd.GetArrayName(i)==fname]
        assert(len(field_numbers)==1)
        if verbose:
            print("output:",pd.GetArrayName(field_numbers[0]) )
        v_interp_on_grid = numpy_support.vtk_to_numpy(pd.GetArray(field_numbers[0]))   
        if shape2d:
            return v_interp_on_grid.reshape(shape2d)
        return v_interp_on_grid

def mod_mesh(du, stlfile='c0006.stl',  output_fn='surface.vtp',
                write=False,sign=1):
    '''
    Move stl in normal direction
    '''
    stl = getFileReaderOutput(stlfile)
    vertices = numpy_support.vtk_to_numpy(stl.GetPoints().GetData())
    indices = numpy_support.vtk_to_numpy(stl.GetPolys().GetData()).reshape(-1, 4)[:, 1:4]
    
    merged = vtk.vtkPolyData()
    merged.DeepCopy(stl)
    

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

    
    points = vtk.vtkPoints()

    for normal, pos in zip(normals, vertices):
        points.InsertNextPoint(pos + normal * (-sign*du))

 
    merged.SetPoints(points)


    return merged

def probe_at_dx(du = 0.001, velocity_file=None, stlfile='c0006.stl',  output_fn='surface.vtp',
                velocity_name = 'v [m/s]',
                write=False, mu=1.0,move_mesh=False,sign=1):
    '''
    Equidistant points from stl in normal direction
    '''
    stl = stlImageActor(stlfile)
    vertices = numpy_support.vtk_to_numpy(stl.GetPoints().GetData())
    indices = numpy_support.vtk_to_numpy(stl.GetPolys().GetData()).reshape(-1, 4)[:, 1:4]
    merged = vtk.vtkPolyData()
    merged.DeepCopy(stl)
    
    vel_data = getFileReaderOutput(velocity_file)

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

    
    points = vtk.vtkPoints()
    pointsPolyData = vtk.vtkPolyData()

    for normal, pos in zip(normals, vertices):
        points.InsertNextPoint(pos + normal * (-sign*du))

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
        

    print(v_vec.shape)

    velocity = vtk.vtkFloatArray()
    velocity.SetNumberOfComponents(1)
    velocity.SetName("X_at_epsilon")

    
    for v in v_vec:

        if np.max(v) > 1e33 or np.min(v) < -1e33:
            s = np.array([ 0.00])

        velocity.InsertNextTypedTuple( [v] )

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

    


def stlImageActor(path):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()

def vtlImageActor(path):
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()

def unstructuredImage(path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    return reader

def getFileReaderOutput(filename):
    import os.path
    _,file_extension = os.path.splitext(filename)
    
    if file_extension.endswith(".vti"):
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(filename)
        
    elif file_extension.endswith(".stl"):
        reader = vtk.vtkSTLReader()
        reader.SetFileName(filename)
        
    elif file_extension.endswith(".vtp"):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(filename)
        
    elif file_extension.endswith(".vtu"):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)

    elif file_extension.endswith(".vtk"):
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        
    reader.Update()
    return reader.GetOutput()




def scale_and_trans(vtk_data=None, output =  None,
                    scale = 1000.0, 
                    deltaxyz = [-14.52326308, 180.637182  , 161.81502267] ):
    '''
    Performs scaling (e.g from meters to mm) and translation of the dataset.
    Note that `vtk_data` is reader.GetOutput()
    '''
    transform = vtk.vtkTransform()
    transform.Scale(scale,scale,scale)
    transform.Translate(*deltaxyz)
    transformFilter = vtk.vtkTransformFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputData(vtk_data)
    transformFilter.Update()
    if output == None:
        return transformFilter.GetOutput()
    else:
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(output)
        writer.SetInputData(transformFilter.GetOutput())
        writer.Update()
        writer.Write()
