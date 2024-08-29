bg_color white
set ray_opaque_background, on
set ambient, 1
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

load PDB/SI06HA_2F01.pdb
create 2F01_HA, SI06HA_2F01
create 2F01_HA1, 2F01_HA and chain C+E+A
create 2F01_HA2, 2F01_HA and chain F+B+D
color grey70, 2F01_HA1
color grey50, 2F01_HA2
create 2F01_H, 2F01_HA and chain H
create 2F01_L, 2F01_HA and chain L
color cb_redorange, 2F01_H
color cb_rose, 2F01_L
distance 2F01_hb_1, 2F01_L and resi 30 and n. OG, 2F01_HA2 and chain B and resi 38 and n. OE1
distance 2F01_hb_2, 2F01_L and resi 32 and n. NE1, 2F01_HA2 and chain B and resi 46 and n. OE1
distance 2F01_hb_3, 2F01_L and resi 94 and n. OH, 2F01_HA2 and chain B and resi 19 and n. OD2
distance 2F01_hb_4, 2F01_L and resi 30 and n. OG, 2F01_HA2 and chain B and resi 38 and n. NE2
color black, 2F01_hb*
set dash_width, 3

hide all
show dash, 2F01_hb*
show cartoon, 2F01_L
show cartoon, 2F01_HA1 and chain A+C or 2F01_HA2 and chain B
show sticks, 2F01_L and resi 30+32+94 and (not name c+n+o) or 2F01_HA2 and chain B and resi 19+38+42 and (not name c+n+o)
util.cnc all
set cartoon_transparency, 0.5, 2F01_HA1
set cartoon_transparency, 0.5, 2F01_HA2
set cartoon_transparency, 0.5, 2F01_L
set_view (\
     0.267632276,   -0.643954873,    0.716715157,\
     0.890518069,   -0.118668780,   -0.439164132,\
     0.367863119,    0.755792618,    0.541697264,\
    -0.007223822,    0.001995094,  -87.906562805,\
   185.655838013,  186.256546021,  144.913909912,\
  -131.683090210,  305.836883545,  -20.000000000 )
ray 1000, 1000
png graph/PDB/2F01_LC.png

load PDB/SI06HA_16ND92.pdb
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
distance 16ND_hb_1, 16ND_L and resi 222 and n. OH, 16ND_HA2 and chain B and resi 38 and n. NE2
distance 16ND_hb_2, 16ND_L and resi 220 and n. ND2, 16ND_HA2 and chain B and resi 42 and n. OE1
distance 16ND_hb_3, 16ND_L and resi 160 and n. NE1, 16ND_HA2 and chain B and resi 46 and n. OD1
distance 16ND_hb_4, 16ND_L and resi 159 and n. OG, 16ND_HA1 and chain D and resi 32 and n. NZ
distance 16ND_hb_5, 16ND_L and resi 222 and n. OH, 16ND_HA2 and chain B and resi 19 and n. O
color gray20, 16ND_hb*
set dash_width, 3

hide all
show cartoon, 16ND_L
show cartoon, 16ND_HA1 and chain A+D or 16ND_HA2 and chain B
show sticks, 16ND_L and resi 159+160+220+222 and (not name c+n+o)
show sticks,  16ND_HA2 and chain B and resi 38+42+46 and (not name c+n+o)
show sticks,  16ND_HA2 and chain B and resi 19
show sticks,  16ND_HA1 and chain D and resi 32 and (not name c+n+o)
show dash, 16ND_hb*
util.cnc all
set cartoon_transparency, 0.5, 16ND_HA1
set cartoon_transparency, 0.5, 16ND_HA2
set cartoon_transparency, 0.5, 16ND_L
set_view (\
     0.267632276,   -0.643954873,    0.716715157,\
     0.890518069,   -0.118668780,   -0.439164132,\
     0.367863119,    0.755792618,    0.541697264,\
    -0.007223822,    0.001995094,  -87.906562805,\
   185.655838013,  186.256546021,  144.913909912,\
  -131.683090210,  305.836883545,  -20.000000000 )
ray 1000, 1000
png graph/PDB/16ND_LC.png
