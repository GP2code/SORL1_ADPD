# Load structure
load q92673_model.cif


# Modify the representation of the system
hide all
preset.publication()

color white, all


# Color domains
select P, resi 29-81
color brown, P

select VPS10, resi 124-757
color salmon, VPS10

select LDLR-B, resi 780-1013
color blue, LDLR-B

select EGF, resi 1026-1074
color orange, EGF

select LDLR-A, resi 1076-1551
color yellow, LDLR-A

select FN3, resi 1555-2112
color pink, FN3

select IC, resi 2159-2214
color purple, IC


# Select the residues of interest
# PD
select K91R    , resi 91    and name CA
select V118M   , resi 118   and name CA
select P171S   , resi 171   and name CA
select F191S   , resi 191   and name CA
select A338V   , resi 338   and name CA
select P345L   , resi 345   and name CA
select N371T   , resi 371   and name CA
select T459K   , resi 459   and name CA
select V718M   , resi 718   and name CA
select A743V   , resi 743   and name CA
select E835D   , resi 835   and name CA
select D864N   , resi 864   and name CA
select R945Q   , resi 945   and name CA
select T947M   , resi 947   and name CA
select R1080C  , resi 1080  and name CA
select F1099L  , resi 1099  and name CA
select F1123L  , resi 1123  and name CA
select R1124S  , resi 1124  and name CA
select Q1157R  , resi 1157  and name CA
select R1172H  , resi 1172  and name CA
select M1279T  , resi 1279  and name CA
select P1454S  , resi 1454  and name CA
select R1490C  , resi 1490  and name CA
select F1519Y  , resi 1519  and name CA
select W1563C  , resi 1563  and name CA
select T1679I  , resi 1679  and name CA
select N1809S  , resi 1809  and name CA
select P1846S  , resi 1846  and name CA
select H1937Q  , resi 1937  and name CA
select L1993F  , resi 1993  and name CA
select D2065V  , resi 2065  and name CA
select R2163Q  , resi 2163  and name CA
select A2184T  , resi 2184  and name CA
select D2190N  , resi 2190  and name CA
select E2194K  , resi 2194  and name CA


# Color the residues of interest
# PD
color magenta, K91R  
color magenta, V118M 
color magenta, P171S 
color magenta, F191S 
color magenta, A338V 
color magenta, P345L 
color magenta, N371T 
color magenta, T459K 
color cyan, V718M
color magenta, A743V 
color magenta, E835D 
color magenta, D864N 
color magenta, R945Q 
color cyan, T947M 
color cyan, R1080C
color cyan, F1099L
color magenta, F1123L
color magenta, R1124S
color cyan, Q1157R
color magenta, R1172H
color magenta, M1279T
color cyan, P1454S
color cyan, R1490C
color magenta, F1519Y
color cyan, W1563C
color magenta, T1679I
color magenta, N1809S
color magenta, P1846S
color magenta, H1937Q
color magenta, L1993F
color cyan, D2065V
color magenta, R2163Q
color magenta, A2184T
color magenta, D2190N
color magenta, E2194K


# Change the representation of the residues of interest
# PD
show spheres, K91R  
show spheres, V118M 
show spheres, P171S 
show spheres, F191S 
show spheres, A338V 
show spheres, P345L 
show spheres, N371T 
show spheres, T459K 
show spheres, V718M 
show spheres, A743V 
show spheres, E835D 
show spheres, D864N 
show spheres, R945Q 
show spheres, T947M 
show spheres, R1080C
show spheres, F1099L
show spheres, F1123L
show spheres, R1124S
show spheres, Q1157R
show spheres, R1172H
show spheres, M1279T
show spheres, P1454S
show spheres, R1490C
show spheres, F1519Y
show spheres, W1563C
show spheres, T1679I
show spheres, N1809S
show spheres, P1846S
show spheres, H1937Q
show spheres, L1993F
show spheres, D2065V
show spheres, R2163Q
show spheres, A2184T
show spheres, D2190N
show spheres, E2194K


# Adjust the size of the label
set label_size=12
set label_position=[7.0, 0.0, 30.0]
set label_font_id=7


# Change the background color
bg_color white
set cartoon_transparency=0.4


# Render image and save
set depth_cue=0
orient
zoom complete=1
#ray 3000, 3000
#png SORL1_PD_labels.png, dpi=600

