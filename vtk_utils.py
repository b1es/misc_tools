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
    # data = data.transpose(2,1,0)
    
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
