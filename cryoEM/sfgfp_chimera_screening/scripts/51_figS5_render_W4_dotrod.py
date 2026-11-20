# Fig S5A-style render, W4_dotrod: RELION loose-mask postprocessed map coloured by RELION local resolution.
# Written by scripts/50_figS5_prep.py; run offscreen:  $CHIMERAX --nogui --offscreen --silent <this file>
from chimerax.core.commands import run
from chimerax.graphics import OrthographicCamera
v = session.main_view
def ortho(fw):
    c = OrthographicCamera(fw)
    c.position = v.camera.position
    v.camera.delete()
    v.camera = c
def frame(fw):
    run(session, 'view #1')               # re-centre on the map after every turn
    v.camera.field_width = fw             # then fix the scale
    v.camera.redraw_needed = True
run(session, 'set bgColor #fcfcfb')
run(session, 'lighting soft')
run(session, 'graphics silhouettes true width 1.5')
run(session, 'open /n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/relion/W4_dotrod/pp_loose/postprocess_masked.mrc')
run(session, 'open /n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/relion/W4_dotrod/relion_locres.mrc')
run(session, 'volume #2 hide')
run(session, 'volume #1 level 0.78705 style surface step 1')
run(session, 'surface dust #1 size 15')
run(session, 'color sample #1 map #2 palette 7.50,#d7191c:8.50,#fdae61:9.50,#ffffbf:10.50,#abdda4:11.50,#2b83ba')
run(session, 'view matrix camera 0.000363,0.079340,-0.996848,-771.087891,0.386985,0.919168,0.073299,286.838164,0.922086,-0.385792,-0.030369,200.369889')
ortho(220.0)
frame(220.0)
run(session, 'save /n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/W4_dotrod/W4_dotrod_v1.png width 1400 height 1100 supersample 3')
run(session, 'turn x 90')
frame(220.0)
run(session, 'save /n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/W4_dotrod/W4_dotrod_v2.png width 1400 height 1100 supersample 3')
run(session, 'turn x -90')
run(session, 'turn y 90')
frame(220.0)
run(session, 'save /n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/W4_dotrod/W4_dotrod_v3.png width 1400 height 1100 supersample 3')
import json
from chimerax.map import Volume
vol = [m for m in session.models.list() if isinstance(m, Volume)][0]
json.dump({'surface_level': float(vol.minimum_surface_level), 'field_width_A': float(v.camera.field_width), 'window_px': [1400, 1100]}, open('/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/W4_dotrod/W4_dotrod_render_log.json', 'w'))
run(session, 'exit')
