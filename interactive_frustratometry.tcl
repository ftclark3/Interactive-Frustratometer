# Establish connection to the persistent Python script
set ::py_pipe [open "|python3 interactive_update.py" r+]
fconfigure $::py_pipe -buffering line -blocking 1 

# Global state trackers
set ::last_physics_pos ""
set ::last_visual_pos ""
set ::cached_frustration_data ""

# =====================================================================
# 1. The Drawing Function
#
#    Whenever any non-glycine CA atom has been moved, 
#    the positions of the lines on the screen may need to be updated, 
#    so we call this function. 
#
#    If a glycine CA atom or any CB atom is moved, then we may need to
#    update not only the positions of the lines on the screen, but also
#    the colors and styles of the lines, possibly adding or deleting 
#    lines. Such adjustments will trigger a recalculation of the 
#    frustratogram, followed by this function.
#
# =====================================================================
proc draw_lines_from_cache {molid} {
    graphics $molid delete all
    
    # If we have no data yet, abort silently
    if {$::cached_frustration_data eq ""} { return }

    # Rebuild CA coordinate array for drawing
    set sel [atomselect $molid "name CA"]
    set atom_data [$sel get {chain resid x y z}]
    $sel delete

    array unset coords
    foreach record $atom_data {
        lassign $record c r x y z
        set coords($c,$r) [list $x $y $z]
    }

    # Parse the cached string from Python
    set lines [split $::cached_frustration_data "\n"]
    
    foreach line $lines {
        set line [string trim $line]
        if {$line eq "" || [string index $line 0] eq "#"} { continue }
        
        set words [split $line " "]
        if {[llength $words] != 14} { continue }
        
        set r1           [lindex $words 0]
        set r2           [lindex $words 1]
        set c1           [lindex $words 2]
        set c2           [lindex $words 3]
        set contact_type [lindex $words 12]
        set state        [lindex $words 13]
        
        if {$state eq "neutral"} { continue }
        
        if {$state eq "minimally"} {
            set color "green"
        } elseif {$state eq "highly"} {
            set color "red"
        } else { continue }

        # solid lines for short-range contacts
        # and dashed lines for long-range contacts 
        # ("water-mediated" and "long," which means "protein-mediated")
        if {$contact_type eq "short"} {
            set style "solid"
        } elseif {$contact_type eq "water-mediated" || $contact_type eq "long"} {
            set style "dashed"
        } else { continue }
        
        if {[info exists coords($c1,$r1)] && [info exists coords($c2,$r2)]} {
            graphics $molid color $color
            graphics $molid line $coords($c1,$r1) $coords($c2,$r2) style $style width 2
        }
    }
}

# =====================================================================
# 2. The Two-Tier Polling Loop
# =====================================================================
proc auto_update_on_move {} {
    set molid 0
    if {[lsearch [molinfo list] $molid] == -1} { return }

    # Get the interaction centers (Physics)
    set sel_phys [atomselect $molid "name CB or (resname GLY IGL and name CA)"]
    set curr_phys [$sel_phys get {x y z}]
    $sel_phys delete

    # Get the visual anchors (Visual)
    set sel_vis [atomselect $molid "name CA"]
    set curr_vis [$sel_vis get {x y z}]
    $sel_vis delete

    set needs_redraw 0

    # TIER 1: Did the physics change?
    if {$curr_phys ne $::last_physics_pos} {
        
        set ::last_physics_pos $curr_phys
        set ::last_visual_pos $curr_vis
        
        # Flatten the coordinate list { {x1 y1 z1} {x2 y2 z2} } -> "x1 y1 z1 x2 y2 z2"
        # and send it straight down the pipe
        puts $::py_pipe [join $curr_phys " "]
        
        # Wait patiently for Python to compute and send data back
        set new_data ""
        while {1} {
            set line [gets $::py_pipe]
            if {$line eq "__END_OF_CALC__"} { break }
            append new_data $line "\n"
        }
        
        # Cache the new data
        set ::cached_frustration_data $new_data
        set needs_redraw 1

    # TIER 2: Did only the atomic positions change?
    } elseif {$curr_vis ne $::last_visual_pos} {
        set ::last_visual_pos $curr_vis
        set needs_redraw 1
    }

    # Redraw only if something actually moved
    if {$needs_redraw} {
        draw_lines_from_cache $molid
    }

    # Reschedule (10ms is fine for updating where the lines are drawn 
    #             but we might need to use a longer period to allow
    #             frustration calculations to finish between line draw updates
    set ::coord_poll_id [after 10 auto_update_on_move]
}

# Clean up and start
catch {after cancel $::coord_poll_id}
auto_update_on_move
