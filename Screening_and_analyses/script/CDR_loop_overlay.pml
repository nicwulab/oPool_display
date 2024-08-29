bg_color white
set ambient, 0.7
set ray_shadows, off
set ray_shadow_fudge, 0
set shininess, 7
set field_of_view, 20
set antialias, 2
set ray_trace_mode, 1
set spec_reflect, 0.35

load PDB/SI06HA_16ND_OO_v2.pdb
create 16ND_mono, SI06HA_16ND_OO_v2 and chain A+B+J+K
hide all
show cartoon,16ND_mono

load PDB/56.a.09_mono.cif
create 56.a.09_mono, chain F or (chain H and resi 1-113) or (chain L and resi 1-108)
align 56.a.09_mono and chain F,  16ND_mono and chain A+B

load PDB/SI06HA2F01v6_real_space_refined.pdb
create 2F01_mono, SI06HA2F01v6_real_space_refined and chain A+B+H+L
align 2F01_mono and chain A+B,  16ND_mono and chain A+B

load PDB/SIA28_mono.cif
align SIA28_mono and chain A+B,  16ND_mono and chain A+B

load PDB/429_B01_mono.cif
create 429_B01, 429_B01_mono and chain A+B or (429_B01_mono and chain H and resi 1-125) or (429_B01_mono and chain L and resi 1-109)
align 429_B01 and chain A+B,  16ND_mono and chain A+B

load PDB/39.29_mono.cif
align 39.29_mono and chain A,  16ND_mono and chain A+B

load PDB/1G05_mono.cif
align 1G05_mono and chain A+B,  16ND_mono and chain A+B

load PDB/MEDI8852_mono.cif
create MEDI8852, MEDI8852_mono and chain A+B or (MEDI8852_mono and chain M and resi 1-127) or (MEDI8852_mono and chain N and resi 1-104)
align MEDI8852 and chain A+B,  16ND_mono and chain A+B

color grey70, 16ND_mono and chain A
color grey50, 16ND_mono and chain B
hide everything

create HA_mono, 16ND_mono and chain A+B

show cartoon, HA_mono
show cartoon, 16ND_mono and chain J and resi 101-111
alter  56.a.09_mono and chain H and resi 95-100B, ss="L"
show cartoon, 56.a.09_mono and chain H and resi 95-100B
alter  2F01_mono and chain H and resi 95-100B, ss="L"
show cartoon, 2F01_mono and chain H and resi  97-100C
alter  SIA28_mono and chain C and resi 130-140, ss="L"
show cartoon, SIA28_mono and chain C and resi 130-138
alter  429_B01 and chain H and resi 99-108, ss="L"
show cartoon, 429_B01 and chain H and resi 99-108
alter  39.29_mono and chain H and resi 100-108, ss="L"
show cartoon, 39.29_mono and chain H and resi 100-108
alter  1G05_mono and chain H and resi 95-100D, ss="L"
show cartoon, 1G05_mono and chain H and resi 96-100B
alter  medi8852 and chain M and resi 103-113, ss="L"
show cartoon, MEDI8852 and chain M and resi 105-111
color marine, MEDI8852
color orange, 16ND_mono

set_view (\
    -0.169601440,   -0.213150576,    0.962180018,\
    -0.943227887,    0.317951500,   -0.095824495,\
    -0.285498202,   -0.923797131,   -0.254973054,\
     0.012384960,    0.013135720,  -50.081851959,\
   203.415649414,  172.401962280,  204.580154419,\
  -3301.805908203, 3403.772216797,  -20.000000000 )

 ray 2000,2000
 png graph/CDR_overlay.png