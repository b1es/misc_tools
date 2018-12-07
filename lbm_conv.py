import vtk
from vtk.util import numpy_support
import numpy as np 


def fix_bbox(cfg):
    if  not 'cuts' in cfg:
        cfg['cuts'] = [[0, 0], [0, 0], [0, 0]]    
        print ("No cuts in config, zeros assumed.")

    x0 = [bb_[0] for bb_ in cfg['bounding_box']]
    padding = np.array(cfg['padding']).reshape(3,2).tolist()
    offset = [c[0]-p[0] for c,p in zip(cfg['cuts'],padding)]

    spacing = [(x[1]-x[0])/(s+sum(c)-sum(p)) for x,s,c,p in zip(cfg['bounding_box'],reversed(cfg['size']),cfg['cuts'],padding)]
    spacing = np.array(spacing)

    # real spacing is rather from longest dim! as given in Voxel size: 0.0473177
    ith = np.argmax(cfg['size'][::-1])
    spacing = np.array([spacing[ith]]*3)

    ### fix box  
    dx = spacing[0]
    cfg['bounding_box'] = [ (x[0]+off*dx +(p[0]-c[0])*dx ,  x[0]+off*dx + (s-(p[0]-c[0]))*dx ) \
     for x,s,c,p,off in zip(cfg['bounding_box'],reversed(cfg['size']),cfg['cuts'],padding,offset)]

    print('cfg fixed!')


def to_phys(pkts,cfg):
    '''
    converts from LB to physical coordinated based on processed config
    uses numpy broadcasting
    attention: LB coordinates are in index order (e.g. as produced from np.where)
    '''
    if  not 'cuts' in cfg:
        cfg['cuts'] = [[0, 0], [0, 0], [0, 0]]    
        print ("No cuts in config, zeros assumed.")
    x0 = [bb_[0] for bb_ in cfg['bounding_box']]
    x_lb = pkts[:,::-1]
    padding = np.array(cfg['padding']).reshape(3,2).tolist()
    offset = [c[0]-p[0] for c,p in zip(cfg['cuts'],padding)]
    spacing = [(x[1]-x[0])/(s+sum(c)-sum(p)) for x,s,c,p in zip(cfg['bounding_box'],reversed(cfg['size']),cfg['cuts'],padding)]
    return (x_lb+offset)*spacing+x0

def get_minmax4voxel(cfg,lattice=False,get_origin_spacing=False):
    '''
    Get en effective/larger bounding box (padding - cutting)
    Parameters:
    
     lattice: default false, if true gives bbox for lattice rather than for voxels.
     get_origin_spacing:  get origin and spacing in cartesian order for e.g. vtk 
    '''
    if lattice:
        lattice = 0.5
    else:
        lattice = 0.0
    if not 'cuts' in cfg.keys():
        print('no cuts in cfg, assuming zeros!')
        cfg['cuts'] = [[0, 0], [0, 0], [0, 0]]    
    padding = np.array(cfg['padding']).reshape(3,2).tolist()
    spacing = [(x[1]-x[0])/(s+sum(c)-sum(p)) for x,s,c,p in zip(cfg['bounding_box'],reversed(cfg['size']),cfg['cuts'],padding)]
    (xmin,xmax),(ymin,ymax),(zmin,zmax) = [(x[0]-dx*(p[0]-c[0]-lattice),x[1]+dx*(p[1]-c[1]-lattice)) for x,c,p,dx in zip(cfg['bounding_box'],cfg['cuts'],padding,spacing)]
    if get_origin_spacing:
        return [xmin,ymin,zmin],spacing
    else:
        return (xmin,xmax),(ymin,ymax),(zmin,zmax)

def prepare_vtk_data(property_name, number_of_components, numpy_data):
    arr = numpy_support.numpy_to_vtk(num_array=numpy_data.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    arr.SetName(property_name)
    arr.SetNumberOfComponents(number_of_components)

    return arr



def convert_to_phys(data,filename="test.vti",cfg=None,drho=1,dx=1,dt=1,visc=-1.0/12.0,verbose=False):
    '''
    Convert to physical unit, geometry (and values)  

    :param cfg: geometry config ( json) 
    :param data: file from sailfish simulation
    '''
    
    if drho==1 and dx==1 and dt==1:
        convert_fields = False
    else:
        convert_fields = True
    
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)

    grid = vtk.vtkImageData()

    if cfg:
        gx,gy,gz = reversed(cfg['size'])
        origin,spacing = get_minmax4voxel(cfg,lattice=True,get_origin_spacing=True)
    else:
        gx,gy,gz = reversed(data['rho'].shape)
        origin = (0,0,0)
        spacing = (1,1,1)
    grid.SetDimensions(gx,gy,gz)
    grid.SetOrigin(*origin)
    grid.SetSpacing(*spacing)

    pd = grid.GetPointData()
    phys_name = dict()
    if hasattr(data,'files'):
        for _name in data.files:
            field = data[_name]
            if convert_fields:
                if _name == 'rho':
                    phys_name[_name] = 'p [Pa]'
                    field *= 1.0/3.0*drho*(dx/dt)**2 
                elif _name == 'v':
                    phys_name[_name] = 'v [m/s]'
                    field *= (dx/dt)
                elif _name.startswith('stress_'):
                    phys_name[_name] =_name + " [Pa]"
                    field *= -6*visc/(1+6*visc)*drho*(dx/dt)**2 
                else:
                    phys_name[_name] =_name
            else:
                phys_name[_name] =_name 

            if verbose:
                print("processing:",_name,phys_name[_name],end='')
            if len(data[_name].shape) == 3:
                if verbose:
                    print(" scalar field:",data[_name].shape,end='')
                pd.AddArray(prepare_vtk_data(phys_name[_name], 1, field))
            elif len(data[_name].shape) == 4:
                assert(data[_name].shape[0]==3)
                field = np.rollaxis(field,0,4)
                pd.AddArray(prepare_vtk_data(phys_name[_name], 3, field)) 
                if verbose:
                    print(" vector field:",data[_name].shape,end='')
            else:
                if verbose:
                    print(" ignoring field with shape:",data[_name].shape,end='')
            if verbose:
                print("... ok")
    else:
        if verbose:
            print("processing a scalar, no units conv, writing to vtk scalar 's' ")
        assert(len(data.shape)==3)
        pd.AddArray(prepare_vtk_data("s", 1, data))
    if verbose:
        print("finished")
    writer.SetInputData(grid)
    writer.Write()
