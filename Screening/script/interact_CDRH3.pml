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
select HP_Groove, 2F01_HA1 and chain A and resi 40+318 or 2F01_HA2 and chain B and resi 21+41+45+48+49+52
alter 2F01_H and chain H and resi 92-102, ss="L"
set cartoon_side_chain_helper, on
hide all
show stick, HP_Groove and (not name c+n+o)
show cartoon, 2F01_HA1 and chain A
show cartoon, 2F01_HA2 and chain B
show cartoon, 2F01_H and chain H and resi 92-102
show sticks,  2F01_H and chain H and resi 100+100A+100B and (not name c+n+o)
set cartoon_transparency, 0.5, 2F01_HA1
set cartoon_transparency, 0.5, 2F01_HA2
set cartoon_transparency, 0.5, 2F01_H
util.cnc all
set_view (\
     0.023823475,    0.031305671,    0.999223053,\
     0.976894438,   -0.213054612,   -0.016615178,\
     0.212369740,    0.976537406,   -0.035658129,\
    -0.001484178,    0.000059120,  -81.887420654,\
   185.281433105,  180.323196411,  152.703186035,\
  -3625.530517578, 3789.385498047,  -20.000000000 )
ray 1000, 1000
png graph/PDB/2F01_CDRH3.png, dpi=300

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
alter 16ND_H and resi 97-115, ss="L"
set cartoon_side_chain_helper, on
select HP_Groove, 16ND_HA1 and chain A and resi 40+42+292+318 or 16ND_HA2 and chain B and resi 21+41+45+48+52+56
hide all
show cartoon, 16ND_HA1 and chain A
show cartoon, 16ND_HA2 and chain B
show sticks, HP_groove and (not name c+n+o)
show cartoon, 16ND_H and resi 97-115
show sticks, 16ND_H and resi 104+105+109+110 and (not name c+n+o)
set cartoon_transparency, 0.5, 16ND_HA1
set cartoon_transparency, 0.5, 16ND_HA2
set cartoon_transparency, 0.5, 16ND_H
util.cnc all
set_view (\
     0.023823475,    0.031305671,    0.999223053,\
     0.976894438,   -0.213054612,   -0.016615178,\
     0.212369740,    0.976537406,   -0.035658129,\
    -0.001577061,   -0.000103924,  -83.078422546,\
   185.262924194,  180.341369629,  152.337921143,\
  -3624.345947266, 3790.570068359,  -20.000000000 )
ray 1000,1000
png graph/PDB/16ND_CDRH3.png, dpi=300
