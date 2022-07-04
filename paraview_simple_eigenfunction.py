# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import os

# create a new 'PVD Reader'
eigenfuction_criticalpvd = PVDReader(registrationName='eigenfuction_critical.pvd', FileName='/Users/fabianlaakmann/Desktop/hartmann/dd/eigenfuction_critical.pvd')
eigenfuction_criticalpvd.PointArrays = ['Velocity', 'Pressure', 'Temperature', 'MagneticField', 'ElectricFieldf']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
eigenfuction_criticalpvdDisplay = Show(eigenfuction_criticalpvd, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'Pressure'
pressureLUT = GetColorTransferFunction('Pressure')

# get opacity transfer function/opacity map for 'Pressure'
pressurePWF = GetOpacityTransferFunction('Pressure')

# trace defaults for the display properties.
eigenfuction_criticalpvdDisplay.Representation = 'Surface'
eigenfuction_criticalpvdDisplay.ColorArrayName = ['POINTS', 'Pressure']
eigenfuction_criticalpvdDisplay.LookupTable = pressureLUT
eigenfuction_criticalpvdDisplay.SelectTCoordArray = 'None'
eigenfuction_criticalpvdDisplay.SelectNormalArray = 'None'
eigenfuction_criticalpvdDisplay.SelectTangentArray = 'None'
eigenfuction_criticalpvdDisplay.OSPRayScaleArray = 'Pressure'
eigenfuction_criticalpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
eigenfuction_criticalpvdDisplay.SelectOrientationVectors = 'Velocity'
eigenfuction_criticalpvdDisplay.ScaleFactor = 0.1
eigenfuction_criticalpvdDisplay.SelectScaleArray = 'Pressure'
eigenfuction_criticalpvdDisplay.GlyphType = 'Arrow'
eigenfuction_criticalpvdDisplay.GlyphTableIndexArray = 'Pressure'
eigenfuction_criticalpvdDisplay.GaussianRadius = 0.005
eigenfuction_criticalpvdDisplay.SetScaleArray = ['POINTS', 'Pressure']
eigenfuction_criticalpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
eigenfuction_criticalpvdDisplay.OpacityArray = ['POINTS', 'Pressure']
eigenfuction_criticalpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
eigenfuction_criticalpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
eigenfuction_criticalpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
eigenfuction_criticalpvdDisplay.ScalarOpacityFunction = pressurePWF
eigenfuction_criticalpvdDisplay.ScalarOpacityUnitDistance = 0.06564197879454707
eigenfuction_criticalpvdDisplay.OpacityArrayName = ['POINTS', 'Pressure']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
eigenfuction_criticalpvdDisplay.ScaleTransferFunction.Points = [-0.02599338207430274, 0.0, 0.5, 0.0, 0.0013370761815960017, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
eigenfuction_criticalpvdDisplay.OpacityTransferFunction.Points = [-0.02599338207430274, 0.0, 0.5, 0.0, 0.0013370761815960017, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

# show color bar/color legend
eigenfuction_criticalpvdDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Rescale transfer function
pressureLUT.RescaleTransferFunction(-0.02599338207430274, 0.001337076181596003)

# Rescale transfer function
pressurePWF.RescaleTransferFunction(-0.02599338207430274, 0.001337076181596003)

# hide color bar/color legend
eigenfuction_criticalpvdDisplay.SetScalarBarVisibility(renderView1, False)

# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

# set scalar coloring
ColorBy(eigenfuction_criticalpvdDisplay, ('POINTS', 'Velocity', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pressureLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
eigenfuction_criticalpvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
eigenfuction_criticalpvdDisplay.SetScalarBarVisibility(renderView1, True)

##
field_dict = {'Velocity': 'u',
              'MagneticField': 'B',
              'Temperature': 'T'}
for i in range(10):
    for field in ['Velocity', 'MagneticField', 'Temperature']:
        # get animation scene
        animationScene1 = GetAnimationScene()

        # Properties modified on animationScene1
        animationScene1.AnimationTime = float(i)

        # get the time-keeper
        timeKeeper1 = GetTimeKeeper()

        # rescale color and/or opacity maps used to exactly fit the current data range
        eigenfuction_criticalpvdDisplay.RescaleTransferFunctionToDataRange(False, True)

        # hide color bar/color legend
        eigenfuction_criticalpvdDisplay.SetScalarBarVisibility(renderView1, False)

        # Properties modified on renderView1
        renderView1.OrientationAxesVisibility = 0

        # set scalar coloring
        ColorBy(eigenfuction_criticalpvdDisplay, ('POINTS', field, 'Magnitude'))

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(pressureLUT, renderView1)

        # rescale color and/or opacity maps used to include current data range
        eigenfuction_criticalpvdDisplay.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        eigenfuction_criticalpvdDisplay.SetScalarBarVisibility(renderView1, True)

        # get color transfer function/color map for field
        velocityLUT = GetColorTransferFunction(field)

        # get opacity transfer function/opacity map for field
        velocityPWF = GetOpacityTransferFunction(field)

        # hide color bar/color legend
        eigenfuction_criticalpvdDisplay.SetScalarBarVisibility(renderView1, False)

        # get layout
        layout1 = GetLayout()

        # layout/tab size in pixels
        layout1.SetSize(1000, 1000)

        # current camera placement for renderView1
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [0.5, 0.5, 10000.0]
        renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
        renderView1.CameraParallelScale = 0.49#0.7071067811865476
        renderView1.ViewSize = [1000, 1000]
        # save screenshot
        save_path = os.path.join(os.getcwd(), f"eigenfunction_{i}_{field_dict[field]}.png")
        SaveScreenshot(save_path, renderView1, ImageResolution=[1000, 1000],
            TransparentBackground=1, 
            # PNG options
            CompressionLevel='0')
