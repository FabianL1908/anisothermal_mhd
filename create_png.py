# trace generated using paraview version 5.10.1
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 10

import os
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

path = os.path.join(os.getcwd(), 'paraview')
branchids = [int(f) for f in os.listdir(path)]

def get_pvd_files(branchid):
    my_path = os.path.join(path, str(branchid))
    files = os.listdir(my_path)
    pvd_files = [f for f in files if f.endswith(".pvd")]
    return pvd_files


for branchid in branchids:
    print(f"Plotting branch {branchid}")
    pvd_files = get_pvd_files(branchid)
    for pvd_file in pvd_files:
        # create a new 'PVD Reader'
        mypvd = PVDReader(registrationName=pvd_file,
                          FileName=os.path.join(path, str(branchid), pvd_file))
        mypvd.PointArrays = ['Velocity', 'Pressure',
                             'Temperature', 'MagneticField', 'ElectricFieldf']

        # get active view
        renderView1 = GetActiveViewOrCreate('RenderView')

        # show data in view
        mypvdDisplay = Show(mypvd, renderView1,
                            'UnstructuredGridRepresentation')

        # get color transfer function/color map for 'Pressure'
        pressureLUT = GetColorTransferFunction('Pressure')

        # get opacity transfer function/opacity map for 'Pressure'
        pressurePWF = GetOpacityTransferFunction('Pressure')

        # trace defaults for the display properties.
        mypvdDisplay.Representation = 'Surface'
        mypvdDisplay.ColorArrayName = ['POINTS', 'Pressure']
        mypvdDisplay.LookupTable = pressureLUT
        mypvdDisplay.SelectTCoordArray = 'None'
        mypvdDisplay.SelectNormalArray = 'None'
        mypvdDisplay.SelectTangentArray = 'None'
        mypvdDisplay.OSPRayScaleArray = 'Pressure'
        mypvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
        mypvdDisplay.SelectOrientationVectors = 'Velocity'
        mypvdDisplay.ScaleFactor = 0.1
        mypvdDisplay.SelectScaleArray = 'Pressure'
        mypvdDisplay.GlyphType = 'Arrow'
        mypvdDisplay.GlyphTableIndexArray = 'Pressure'
        mypvdDisplay.GaussianRadius = 0.005
        mypvdDisplay.SetScaleArray = ['POINTS', 'Pressure']
        mypvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
        mypvdDisplay.OpacityArray = ['POINTS', 'Pressure']
        mypvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
        mypvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
        mypvdDisplay.PolarAxes = 'PolarAxesRepresentation'
        mypvdDisplay.ScalarOpacityFunction = pressurePWF
        mypvdDisplay.ScalarOpacityUnitDistance = 0.06564197879454707
        mypvdDisplay.OpacityArrayName = ['POINTS', 'Pressure']

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        mypvdDisplay.ScaleTransferFunction.Points = [
            -34174.37950667546, 0.0, 0.5, 0.0, 15990.774417215282, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        mypvdDisplay.OpacityTransferFunction.Points = [
            -34174.37950667546, 0.0, 0.5, 0.0, 15990.774417215282, 1.0, 0.5, 0.0]

        # reset view to fit data
        renderView1.ResetCamera(False)

        # changing interaction mode based on data extents
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [0.5, 0.5, 10000.0]
        renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

        # show color bar/color legend
        mypvdDisplay.SetScalarBarVisibility(renderView1, True)

        # update the view to ensure updated data information
        renderView1.Update()

        # Rescale transfer function
        pressureLUT.RescaleTransferFunction(-34174.37950667546,
                                            15990.774417215285)

        # Rescale transfer function
        pressurePWF.RescaleTransferFunction(-34174.37950667546,
                                            15990.774417215285)

        # Properties modified on renderView1
        renderView1.OrientationAxesVisibility = 0

        # set scalar coloring
        ColorBy(mypvdDisplay, ('POINTS', 'Velocity', 'Magnitude'))

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(pressureLUT, renderView1)

        # rescale color and/or opacity maps used to include current data range
        mypvdDisplay.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        mypvdDisplay.SetScalarBarVisibility(renderView1, True)

        # get color transfer function/color map for 'Velocity'
        velocityLUT = GetColorTransferFunction('Velocity')

        # get opacity transfer function/opacity map for 'Velocity'
        velocityPWF = GetOpacityTransferFunction('Velocity')

        # hide color bar/color legend
        mypvdDisplay.SetScalarBarVisibility(renderView1, False)

        # get layout
        layout1 = GetLayout()

        # layout/tab size in pixels
        layout1.SetSize(1000, 1000)

        # current camera placement for renderView1
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [0.5, 0.5, 10000.0]
        renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
        renderView1.CameraParallelScale = 0.49

        renderView1.ViewSize = [1000, 1000]
        # save screenshot
        SaveScreenshot(os.path.join(path, str(branchid), pvd_file.split(".")[0]+"_u.png"), renderView1, ImageResolution=[1000, 1000],
                       TransparentBackground=1, CompressionLevel='3')

        # set scalar coloring
        ColorBy(mypvdDisplay, ('POINTS', 'Temperature'))

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(velocityLUT, renderView1)

        # rescale color and/or opacity maps used to include current data range
        mypvdDisplay.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        mypvdDisplay.SetScalarBarVisibility(renderView1, True)

        # get color transfer function/color map for 'Temperature'
        temperatureLUT = GetColorTransferFunction('Temperature')

        # get opacity transfer function/opacity map for 'Temperature'
        temperaturePWF = GetOpacityTransferFunction('Temperature')

        # hide color bar/color legend
        mypvdDisplay.SetScalarBarVisibility(renderView1, False)

        # layout/tab size in pixels
        layout1.SetSize(1000, 1000)

        # current camera placement for renderView1
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [0.5, 0.5, 10000.0]
        renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
        renderView1.CameraParallelScale = 0.5

        # save screenshot
        SaveScreenshot(os.path.join(path, str(branchid), pvd_file.split(".")[0]+"_T.png"), renderView1, ImageResolution=[1000, 1000],
                       TransparentBackground=1,
                       # PNG options
                       CompressionLevel='3')

        # set scalar coloring
        ColorBy(mypvdDisplay, ('POINTS', 'MagneticField', 'Magnitude'))

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(temperatureLUT, renderView1)

        # rescale color and/or opacity maps used to include current data range
        mypvdDisplay.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        mypvdDisplay.SetScalarBarVisibility(renderView1, True)

        # get color transfer function/color map for 'MagneticField'
        magneticFieldLUT = GetColorTransferFunction('MagneticField')

        # get opacity transfer function/opacity map for 'MagneticField'
        magneticFieldPWF = GetOpacityTransferFunction('MagneticField')

        # hide color bar/color legend
        mypvdDisplay.SetScalarBarVisibility(renderView1, False)

        # layout/tab size in pixels
        layout1.SetSize(1516, 968)

        # current camera placement for renderView1
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [0.5, 0.5, 10000.0]
        renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
        renderView1.CameraParallelScale = 0.5

        # save screenshot
        SaveScreenshot(os.path.join(path, str(branchid), pvd_file.split(".")[0]+"_B.png"), renderView1, ImageResolution=[1000, 1000],
                       TransparentBackground=1,
                       # PNG options
                       CompressionLevel='3')

        # ================================================================
        # addendum: following script captures some of the application
        # state to faithfully reproduce the visualization during playback
        # ================================================================

        # --------------------------------
        # saving layout sizes for layouts

        # layout/tab size in pixels
        layout1.SetSize(1000, 1000)

        # -----------------------------------
        # saving camera placements for views

        # current camera placement for renderView1
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [0.5, 0.5, 10000.0]
        renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
        renderView1.CameraParallelScale = 0.5

        # --------------------------------------------
        # uncomment the following to render all views
        # RenderAllViews()
        # alternatively, if you want to write images, you can use SaveScreenshot(...).
