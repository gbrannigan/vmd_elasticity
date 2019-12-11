;# time notes here are built by a much larger system than the bench marks
;# all times were arrived at by multple trials with a given number of frames and beads
proc vecavg {list_in} {
    return [expr 1.0*[vecsum $list_in]/[llength $list_in]]
}
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
proc main {lipidinput hg end} {
    ;# takes ~ 3.5 seconds to go through 1 resid and 1000 frames
    ;# takes ~ 25 seconds to run through ~4000 lipids and a 1 frame (~6 sec)
    ;# a 40x40 nm system ~ 5000 lipids will take ~ 5-6 hours to analyze
    ;# MAIN
    set res_sel [atomselect top "$lipidinput"] 
    set resid_list [lsort -unique [$res_sel get resid]]
    $res_sel delete
    set f [open "/u1/home/lms464/lms464/github/JPC_Special/tasks/10_Monge/data/test_output.dat" w]
    set nframes [molinfo top get numframes]
    foreach rid $resid_list {
        set res_sel [atomselect top "resid ${rid} and name $hg" ]
        set temp_sel [atomselect top "(name $hg) and (pbwithin 50 of (resid ${rid} and name $hg))" ]
        set resids_leaflet_list {}
        for {set frm 0} {$frm < $end} {incr frm 1} {
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
