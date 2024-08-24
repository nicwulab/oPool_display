bg_color white
set ray_opaque_background, on
set ambient, 0.7
set ray_shadows, off
set ray_shadow_fudge, 0
set shininess, 7
set field_of_view, 20
set antialias, 2
set ray_trace_mode, 1
set spec_reflect, 0.35
set valence, 0
bg_color white
run https://github.com/Pymol-Scripts/Pymol-script-repo/raw/master/colorblindfriendly.py

load data/PDB/SI06HA_2F01.pdb
create 2F01_HA, SI06HA_2F01
create 2F01_HA1, 2F01_HA and chain C+E+A
create 2F01_HA2, 2F01_HA and chain F+B+D
color grey70, 2F01_HA1
color grey50, 2F01_HA2
create 2F01_H, 2F01_HA and chain H
create 2F01_L, 2F01_HA and chain L
color cb_redorange, 2F01_H
color cb_rose, 2F01_L
hide all
show surface, 2F01_HA1
show surface, 2F01_HA2
show cartoon, 2F01_H
show cartoon, 2F01_L
set_view (\
     0.376762390,   -0.031921051,    0.925757468,\
     0.926219463,   -0.000838214,   -0.376980782,\
     0.012810591,    0.999488652,    0.029250216,\
    -0.000121534,    0.000430785, -482.604614258,\
   174.077667236,  177.757766724,  165.414382935,\
   302.987976074,  662.217468262,  -20.000000000 )
ray 1000, 1000
png graph/PDB/2F01_overview.png

load data/PDB/SI06HA_16ND92.pdb
create 16ND_HA, SI06HA_16ND92
align 16ND_HA and chain A+B+C+D+E+F, 2F01_HA and chain A+B+C+D+E+F
create 16ND_HA1, 16ND_HA and chain A+D+E
create 16ND_HA2, 16ND_HA and chain B+C+F
color grey70, 16ND_HA1
color grey50, 16ND_HA2
create 16ND_H, 16ND_HA and chain J
create 16ND_L, 16ND_HA and chain K
color cb_redorange, 16ND_H
color cb_rose, 16ND_L
hide all
show surface, 16ND_HA1
show surface, 16ND_HA2
show cartoon, 16ND_H
show cartoon, 16ND_L
set_view (\
     0.376762390,   -0.031921051,    0.925757468,\
     0.926219463,   -0.000838214,   -0.376980782,\
     0.012810591,    0.999488652,    0.029250216,\
    -0.000121534,    0.000430785, -482.604614258,\
   174.077667236,  177.757766724,  165.414382935,\
   302.987976074,  662.217468262,  -20.000000000 )
ray 1000,1000
png graph/PDB/16ND_overview.png

