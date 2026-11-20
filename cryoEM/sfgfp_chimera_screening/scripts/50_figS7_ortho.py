# Switch the ChimeraX camera to orthographic (offscreen 'camera ortho' fails on OffScreenRenderingContext.stereo).
from chimerax import graphics
v = session.main_view  # noqa: F821  (session is provided by ChimeraX 'open')
v.camera = graphics.OrthographicCamera(v.camera.view_width(v.center_of_rotation))
