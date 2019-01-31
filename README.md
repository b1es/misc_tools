## Tips and tricks

#### Adding this repository to a project

The easiest way is to use git's submodules feature:
```git
git submodule add git@github.com:marcinofulus/misc_tools.git misc_tools
```

#### Logging vs print and jupyter notebooks

Because this module can be imported to some non-jupyter project, it's good 
practice to use the `logging` module instead of raw `print()` calls.

To make log entries visible in jupyter notebooks, copy the following code:
```python
import logging
logging.getLogger().setLevel(logging.DEBUG)
```

To add timestamps (with milliseconds resolution) to those entries, use:
```python
import logging
formatter = logging.Formatter('[%(asctime)s] [%(levelname)7s] [%(name)s] %(message)s')
for handler in logging.getLogger().handlers:
    handler.setFormatter(formatter)
```

## Resources

* VTK Textbook  
  https://vtk.org/vtk-textbook/  
  https://gitlab.kitware.com/vtk/textbook
* VTK User's Guide  
  https://vtk.org/vtk-users-guide/


## Changelog

#### 2019.01.31 
Renaming:
* `vtk2npy` --> `sailfish_vti_to_npy`
* `stlImageActor` --> `read_vtk`
* `vtlImageActor` --> `read_vtk`
* `unstructuredImage` --> `read_vtk`
* `getFileReaderOutput` --> `read_vtk`

Logging / print
* replaced all `print()` calls with `logging` module.

Added:
* conversions to and from `numpy`
* `SimpleITK` conversions