;# Liam Sharp Dec 2019
;# 
;#		READ ME
;# 1) Assumption VMD is open with no system loaded
;#      -call load_system "structure file" "trajectory file"
;#	  this will load both files for analysis
;#      -call main "lipid species" "head group"
;#	  lipid species tells the function which lipid to look at
;#	  head group tells the function where to measure the 
;#	  leaflet surface
;#
;#		EXAMPLES
;#
;# 2) Example Coarse Grained, assumption I am working on a linux or mac
;# 	(base)$ vmd
;#	% load_system "memb_b_01_wrapped.gro" "memb_b_01_wrapped.xtc"
;#	% main "resname DUPE" "PO4"
;#
;#
;# 3) Example Atomistic (tested but not nessisarily stable)
;#	(base)$ vmd
;#	% load_system "memb_b_01_wrapped.gro" "memb_b_01_wrapped.xtc"
;#	% main "resname DUPE" "P" ;# phosphorous of the PO4
;#
;#
;#		Extra Notes
;#
;# I have a number of other functions that have been tested and have 
;#  consistently worked (at least on my work station)
;#	-Center_System ;# centers the system
;#	-Get_Last_Acyl_Bead ;# determines what the final acyl chain beads are
;#	-no_lipids ;# test to confirm lipids are in the system
;#	-aa_cg_comp ;# attempts to determine if you are using
;#		atomistic or coarse grained (not 100% stable)
;#
;#
;# 		Functions
;#
;# time notes here are built by a much larger system than the bench marks
;# all times were arrived at by multple trials with a given number of frames and beads

proc load_system {struct traj} {
     ;# left with no file type incase user uses gro, pdb, xtc, dcd...	    
     mol new "${struct}"
     mol addfile "${traj}" first 0 last -1 step 1 waitfor -1
}

;# returns a scarlar sum of an the list input
proc vecavg {list_in} {
    return [expr 1.0*[vecsum $list_in]/[llength $list_in]]
}

;# Determines if a lipid is in the extracellular 
;#   or intercellular leaflet 

proc local_mid_plane {atsel_in temp_in} {
    ;# get z ~ 220 us, get z with vecavg ~315 us
    ;# measure center and lindex z ~ 340 us
    ;# saves ~ 25 us using current method
    set resid_z [${atsel_in} get z]; set mid_point [vecavg [${temp_in} get z] ]    
    if {$mid_point < $resid_z} {
        return 0
    } else {
        return 1
    }
}
;# Main and local_mid_plane are the two functions called here
proc main {lipidinput hg} {
    ;# takes ~ 3.5 seconds to go through 1 resid and 1000 frames
    ;# takes ~ 25 seconds to run through ~4000 lipids and a 1 frame (~6 sec)
    ;# a 40x40 nm system ~ 5000 lipids will take ~ 5-6 hours to analyze
    ;# MAIN
    set res_sel [atomselect top "$lipidinput"] 
    set resid_list [lsort -unique [$res_sel get resid]]
    $res_sel delete
    set f [open "test_output.dat" w]
    set nframes [molinfo top get numframes]
    foreach rid $resid_list {
	;# reference lipid
        set res_sel [atomselect top "resid ${rid} and name $hg" ]
	;# surounding lipids, used in local_mid_plane
        set temp_sel [atomselect top "(name $hg) and (pbwithin 50 of (resid ${rid} and name $hg))" ]
        set resids_leaflet_list {}
        for {set frm 0} {$frm < $nframes} {incr frm 1} {
	    ;# update and move frames forwards
            ;# update 50 us faster than atomselect and delete..
            $res_sel frame ${frm}; $temp_sel frame ${frm}; $temp_sel update
            ;# instead of constantly saving data to an array and as set of lists
            ;# ~170 us to run local_mid_plane2 and append, a few us slower to do it one line
            set tmps [local_mid_plane ${res_sel} ${temp_sel}]
            lappend resids_leaflet_list $tmps
        }
        $res_sel delete
        $temp_sel delete
        puts $f "${rid}\t${resids_leaflet_list}"
    }
    close $f
}
