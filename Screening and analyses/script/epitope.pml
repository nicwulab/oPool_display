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

load data/PDB/SI06HA_2F01.pdb
create 2F01_HA, SI06HA_2F01
create 2F01_HA1, 2F01_HA and chain C+E+A
create 2F01_HA2, 2F01_HA and chain F+B+D
color grey70, 2F01_HA1
color grey50, 2F01_HA2
select VH_epitope, (2F01_HA1 and chain A and resi 18+20+37+38+40+318+319 or 2F01_HA2 and chain B and resi 16-18+20,21,41,48)
color cb_red_orange, VH_epitope
select VL_epitope, (2F01_HA2 and chain B and resi 36+39+43+46+49+150+154)
color cb_rose, VL_epitope
select shared_epitope, (2F01_HA2 and chain B and resi 19+34+35+38+42+45)
color violetpurple, shared_epitope
hide all
show surface, 2F01_HA1
show surface, 2F01_HA2
set_view (\
-0.483078837,   -0.112589620,    0.868303478,\
0.875570714,   -0.064122640,    0.478807032,\
0.001768433,    0.991561830,    0.129557610,\
-0.000724005,   -0.001854211, -262.151062012,\
168.970535278,  170.107391357,  138.100692749,\
-2433.562011719, 2957.605224609,  -20.000000000 )
ray 1000,1000
png graph/PDB/2F01_epitope.png

load data/PDB/SI06HA_16ND92.pdb
create 16ND_HA, SI06HA_16ND92
align 16ND_HA and chain A+B+C+D+E+F, 2F01_HA and chain A+B+C+D+E+F
create 16ND_HA1, 16ND_HA and chain A+D+E
create 16ND_HA2, 16ND_HA and chain B+C+F
color grey70, 16ND_HA1
color grey50, 16ND_HA2
select VH_epitope, (16ND_HA1 and chain A and resi 18+38+39+40+42+291-293+318 or 16ND_HA2 and chain B and resi 18-21+41+45+48+49+52+53+56)
color cb_redorange, VH_epitope
select VL_epitope, (16ND_HA1 and chain A and resi 280 or 16ND_HA1 and chain D and resi 32+33 or 16ND_HA2 and chain B and resi 19+20+38+41+42+45+46+49+50+53+57)
color cb_rose, VL_epitope
select shared_epitope, (16ND_HA2 and chain B and resi 19+20+41+45+49)
color violetpurple, shared_epitope
hide all
show surface, 16ND_HA1
show surface, 16ND_HA2
set_view (\
-0.483078837,   -0.112589620,    0.868303478,\
0.875570714,   -0.064122640,    0.478807032,\
0.001768433,    0.991561830,    0.129557610,\
-0.000724005,   -0.001854211, -262.151062012,\
168.970535278,  170.107391357,  138.100692749,\
-2433.562011719, 2957.605224609,  -20.000000000 )
ray 1000,1000
png graph/PDB/16ND_epitope.png

