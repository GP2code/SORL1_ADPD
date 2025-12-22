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
# ADRD
select M1T     , resi 1     and name CA
select M105T   , resi 105   and name CA
select R303W   , resi 303   and name CA
select Y391C   , resi 391   and name CA
select R416Q   , resi 416   and name CA
select V623A   , resi 623   and name CA
select P885S   , resi 885   and name CA
select R985X   , resi 985   and name CA
select C1030fs , resi 1030  and name CA
select R1084C  , resi 1084  and name CA
select D1108N  , resi 1108  and name CA
select R1286C  , resi 1286  and name CA
select D1348G  , resi 1348  and name CA
select C1368R  , resi 1368  and name CA
select R1866W  , resi 1866  and name CA
select F1873Y  , resi 1873  and name CA
select G2090V  , resi 2090  and name CA
select A2171T  , resi 2171  and name CA


# Color the residues of interest
# ADRD
color magenta, M1T    
color cyan, M105T  
color cyan, R303W  
color cyan, Y391C  
color magenta, R416Q  
color magenta, V623A  
color magenta, P885S  
color magenta, R985X  
color magenta, C1030fs
color cyan, R1084C 
color magenta, D1108N 
color magenta, R1286C 
color magenta, D1348G 
color magenta, C1368R 
color magenta, R1866W 
color magenta, F1873Y 
color cyan, G2090V 
color magenta, A2171T 


# Change the representation of the residues of interest
# ADRD
show spheres, M1T    
show spheres, M105T  
show spheres, R303W  
show spheres, Y391C  
show spheres, R416Q  
show spheres, V623A  
show spheres, P885S  
show spheres, R985X  
show spheres, C1030fs
show spheres, R1084C 
show spheres, D1108N 
show spheres, R1286C 
show spheres, D1348G 
show spheres, C1368R 
show spheres, R1866W 
show spheres, F1873Y 
show spheres, G2090V 
show spheres, A2171T 


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
#png SORL1_ADRD_labels.png, dpi=600

