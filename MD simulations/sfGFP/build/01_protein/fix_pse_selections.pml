# Corrected attachment-site selections for the sfGFP-DNA spring chimera.
# Load over the existing session:  pymol "sfGFP 2b3p (150TAG, 134TAG).pse" fix_pse_selections.pml
#
# The session's 134TAG / 150TAG selections point at 2B3P Gly134 and Val150.
# Those are one residue off from the sites actually encoded.  The construct is
# numbered +1 relative to the 2B3P deposition, and the ESI-MS shifts settle it:
#   single Tet construct  27,968 - 27,827 = +141 Da  = Asn -> Tet2-Et
#   double Tet construct  28,108 - 27,827 = +281 Da  = +141 (Asn) + 140 (Asp)
# Gly -> Tet2-Et would be +198 Da and Val -> Tet2-Et +156 Da; neither is observed.
delete 134TAG_wrong
delete 150TAG_wrong
select 134TAG_wrong, 2b3p and chain A and resi 134     # Gly134, the old bookmark
select 150TAG_wrong, 2b3p and chain A and resi 150     # Val150, the old bookmark
select site_D134, 2b3p and chain A and resi 133        # Asp133 in 2B3P = construct D134
select site_N150, 2b3p and chain A and resi 149        # Asn149 in 2B3P = construct N150
select His148_gate, 2b3p and chain A and resi 148      # donates to the chromophore phenolate
show sticks, site_D134 or site_N150 or His148_gate or chromo
color orange, site_D134
color marine, site_N150
color yellow, His148_gate
distance span, site_D134 and name CB, site_N150 and name CB
print "Asp133(CB)-Asn149(CB) = 31.3 A"
