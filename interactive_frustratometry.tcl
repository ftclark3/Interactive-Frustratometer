# SCRIPT TO CALCULATE FRUSTRATION ON-THE-FLY AS COORDINATES
# OF A STRUCTURE ARE CHANGED.
#
# IF MULTIPLE MOLECULES ARE LOADED,
# THE FRUSTRATOGRAM WILL BE CALCULATED FOR THE TOP MOLECULE
# AND RE-CALCULATED AS THE TOP MOLECULE IS CHANGED.
# THIS REQUIRES CLOSING AND RE-LAUNCHING THE PYTHON PIPE,
# WHICH WILL CAUSE A DELAY OF UP TO A FEW (~10) SECONDS.
#
# ONCE THE PYTHON PIPE IS INITIALIZED, RECALCULATION OF THE
# FRUSTRATOGRAM FOLLOWING AN ADJUSTMENT OF MOLECULAR COORDINATES
# IS ORDERS OF MAGNITUDE FASTER.
#
# YOU MAY BE PROMPTED FOR THE SEQUENCE AND LOCATIONS OF THE
# CHAIN BREAKS IN YOUR STRUCTURE IF THE SEQUENCE CANNOT BE
# DETERMINED FROM YOUR STRUCTURE FILE (E.G., IN THE CASE
# THAT YOUR STRUCTURE FILE USES IGL/IPR/NGP RESIDUE NAMES
# INSTEAD OF STANDARD 3-LETTER RESIDUE NAMES)


# =====================================================================
# Safe Tkinter Pop-up for User Input
# (need to ask the user for the sequence and chain breaks if
#  first attempt at structure initialization fails)
# =====================================================================
proc ask_user_for_python {prompt_text} {
    package require Tk
    set ::user_response ""

    set t [toplevel .py_prompt]
    wm title $t "Frustratometer: Input Required"
    wm attributes $t -topmost 1 ;# Force window above VMD

    label $t.lbl -text $prompt_text -wraplength 400 -justify left
    entry $t.ent -textvariable ::user_response -width 60
    button $t.btn -text "Submit" -command {destroy .py_prompt}

    pack $t.lbl $t.ent $t.btn -padx 15 -pady 15

    # Let the user press Enter to submit
    focus $t.ent
    bind $t.ent <Return> {destroy .py_prompt}

    # tkwait pauses the script until the window is destroyed,
    # but keeps the VMD graphics loop alive!
    tkwait window $t
    return $::user_response
}


# =====================================================================
# Procedure to terminate this script,
# called in a few different cases
# (CTRL-F kill_frustration_script for details)
# =====================================================================
proc kill_frustration_script {reason} {
    puts "\n=========================================================="
    puts "FRUSTRATOMETER SHUTDOWN: $reason"
    puts "==========================================================\n"

    # Close the pipe if it exists
    if {[info exists ::py_pipe]} {
        catch {close $::py_pipe}
    }

    # Clear any drawn lines on the active molecule
    if {$::active_molid != -1 && [lsearch [molinfo list] $::active_molid] != -1} {
        catch {graphics $::active_molid delete all}
    }

    # This procedure is called from the polling loop when it is determined that the
    # user has altered the topology. The actual termination is handled
    # in the calling scope, which throws an early return, causing the polling loop
    # to end and this script to die naturally.
}


# =====================================================================
# The Context Switcher (Reads from disk)
# =====================================================================
proc start_python_for_mol {molid} {
    if {[info exists ::py_pipe]} { catch {close $::py_pipe} }

    # Extract the original file path.
    # 'molinfo get filename' returns a nested list: {{file1.pdb} {file2.dcd}}
    # We grab the very first file, which is almost always the structure file.
    set struct_file [lindex [lindex [molinfo $molid get filename] 0] 0]

    # Verify the file still exists on disk
    if {![file exists $struct_file]} {
        kill_frustration_script "Cannot locate original structure file on disk: $struct_file"
        return 0 ;# Return failure flag
    }

    puts "Starting Python server for molecule $molid using file: $struct_file"
    # helpful exception handling in the event that python fails to launch
    if {[catch {open "|python3 -u update_frustration.py \"$struct_file\"" r+} ::py_pipe]} {
        kill_frustration_script "OS level failure: $::py_pipe"
        return 0
    }

    fconfigure $::py_pipe -buffering line -blocking 1

    # We wait here until Python says it is ready, or asks for user input
    # (which is collected using our ask_user_for_python procedure)
    while {1} {
        set char_count [gets $::py_pipe line]

        # Check if the channel hit EOF (Python crashed during initialization)
        if {[eof $::py_pipe] || $char_count == -1} {
            # Close the pipe inside a catch block to extract the Python traceback
            catch {close $::py_pipe} python_error
            kill_frustration_script "Python failed to start.\nTraceback/Error:\n$python_error"
            return 0
        }

        if {$line eq "__READY__"} {
            puts "Frustratometer initialized and ready!"
            break
        } elseif {$line eq "__PROMPT__"} {
            # Python wants input. The very next line in the pipe is the question.
            gets $::py_pipe prompt_text

            # Show the pop-up and wait for user to hit submit
            set user_val [ask_user_for_python $prompt_text]

            # Send the answer back down the pipe
            puts $::py_pipe $user_val
            flush $::py_pipe
        } elseif {[string match "#*" $line]} {
            # Print any debug statements Python sends during init
            puts $line
        }
    }

    set ::active_molid $molid
    set ::last_physics_pos ""
    set ::last_visual_pos ""
    set ::cached_frustration_data ""

    return 1 ;# Return success flag
}


# =====================================================================
#    The Drawing Function
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


    # clear any previously drawn lines
    graphics $molid delete all


    # If we have no data yet, abort silently
    if {$::cached_frustration_data eq ""} { return }


    # get information about the lines that need to be drawn
    set lines [split $::cached_frustration_data "\n"]


    # build index of all CA coordinates using (chain, resid) pairs
    # (these coordinates are needed to determine the endpoints of the lines)
    set sel [atomselect $molid "name CA"]
    #set atom_data [$sel get {residue x y z}]
    set atom_data [$sel get {chain resid x y z}]
    $sel delete
    array unset coords
    foreach record $atom_data {
        lassign $record chain resid x y z
        set chain [string trim $chain]
        set coords($chain,$resid) [list $x $y $z]
    }


    # loop over per-line information, drawing lines
    # based on the CA coordinates and frustration data that we just calculated
    foreach line $lines {

        # if python caught an exception and returned an error message
        # instead of the normal data, print the error and exit
        if {[string match "# PYTHON ERROR*" $line]} {
            puts "\n=== ERROR CAUGHT IN PYTHON ==="
            puts $::cached_frustration_data
            puts "==============================\n"
            continue
        }

        # ignore any comment lines, which we usually don't expect to be there
        if {[string index $line 0] eq "#"} { continue }

        # gets rid of trailing newline?
        set line [string trim $line]

        # ignore any blank lines, which we usually don't expect to be there
        if {$line eq ""} { continue }

        # assume that we have a normal line of data formatted like:
        #     f'{Res1s[counter]} {Res2s[counter]} A A foo foo foo foo foo foo foo foo {style} {descriptor}\n'
        set words [split $line " "]
        #if {[llength $words] != 14} { continue }

        # Grab Chain and Residue IDs from columns 0-3 and contact type and frustration level from 12-13
        set resid1   [lindex $words 0]
        set resid2   [lindex $words 1]
        set c1       [string trim [lindex $words 2]]
        set c2       [string trim [lindex $words 3]]
        set contact_type [lindex $words 12]
        set state        [lindex $words 13]

        # determine whether/how to draw line
        if {$state eq "neutral"} { continue }
        if {$state eq "minimally"} {
            set color "green"
        } elseif {$state eq "highly"} {
            set color "red"
        } else { continue }
        if {$contact_type eq "short" || $contact_type eq "long"} {
            set style "solid"
        } elseif {$contact_type eq "water-mediated"} {
            set style "dashed"
        } else { continue }

        # draw the line
        if {[info exists coords($c1,$resid1)] && [info exists coords($c2,$resid2)]} {
            graphics $molid color $color
            #graphics $molid line $coords($res1_idx) $coords($res2_idx) style $style width 2
            graphics $molid line $coords($c1,$resid1) $coords($c2,$resid2) style $style width 2
        }
    }
}

# =====================================================================
#  The Two-Tier Polling Loop:
#  Periodically checks whether coordinates have changed and
#  moves the locations of the ends of the lines drawn on the screen,
#  first recalculating which lines need to be drawn, if necessary.
# =====================================================================
proc auto_update_on_move {} {

    # -----------------------------------------------------------------
    # THE TOPOLOGY WATCHER (Kill Switch)
    # -----------------------------------------------------------------
    foreach m [molinfo list] {
        set current_atoms [molinfo $m get numatoms]

        # Record atom count if we haven't seen this molecule yet
        if {![info exists ::mol_initial_atoms($m)]} {
            set ::mol_initial_atoms($m) $current_atoms
        }

        # Trip the kill switch if it changed!
        # I guess this script will fail to catch cases where we change
        # an atom name/element but that's far outside of what I would
        # expect the user to try to do.
        if {$current_atoms != $::mol_initial_atoms($m)} {
            kill_frustration_script "Topology changed (atom count mismatch) in molecule $m.\n We are not sure that we can compute frustration correctly\n when the topology of a VMD entry is modified within a session. Out of an abundance of caution, this interactive frustration handler tcl script will be terminated. If you want to compare two similar structures, you can load them in VMD as different entries."
            return ;# IMMEDIATELY exit so the loop cannot reschedule itself
        }
    }

    # -----------------------------------------------------------------
    # THE CONTEXT SWITCHER (resets python if we change our molecule)
    # -----------------------------------------------------------------
    set top_mol [molinfo top]

    if {$top_mol == -1} {
        set ::coord_poll_id [after 10 auto_update_on_move]
        return
    }

    # If the user switched to a different molecule
    if {$top_mol != $::active_molid} {

        # Clear lines on the old molecule before switching
        if {$::active_molid != -1 && [lsearch [molinfo list] $::active_molid] != -1} {
            catch {graphics $::active_molid delete all}
        }

        # If start_python_for_mol returns 0 (fails), abort the loop entirely
        if {![start_python_for_mol $top_mol]} {
            return
        }

    }

    set molid $::active_molid


    # -----------------------------------------------------------------
    # CHECK IF WE NEED TO REDRAW OUR LINES, WITH OR WITHOUT
    # A RECALCULATION OF THE FRUSTRATOGRAM.
    # TYPICALLY, IT WILL BE NECESSARY TO RECALCULATE THE FRUSTRATOGRAM,
    # AND WE HAVE OPTIMIZED THE SPEED OF THIS CALCULATION FOR THE
    # PURPOSES OF THIS INTERACTIVE SCRIPT.
    # ----------------------------------------------------------------
    # Get current state
    set current_frame [molinfo $molid get frame]
    set struct_file [lindex [lindex [molinfo $molid get filename] 0] 0]
    
    # Create a safe, unique filename (e.g., myprotein.pdb_frame_5.dat)
    set safe_name [file tail $struct_file]
    set cache_file [file join $::frust_cache_dir "${safe_name}_frame_${current_frame}.dat"]

    # Get the interaction centers and visual anchors
    set sel_phys [atomselect $molid "name CB or (resname GLY IGL and name CA)"]
    set curr_phys [$sel_phys get {x y z}]
    $sel_phys delete

    set sel_vis [atomselect $molid "name CA"]
    set curr_vis [$sel_vis get {x y z}]
    $sel_vis delete

    set needs_redraw 0

    # -----------------------------------------------------------------
    # CACHE CHECK: Did the frame change?
    # -----------------------------------------------------------------
    if {$current_frame != $::last_frame} {
        set ::last_frame $current_frame
        
        if {[file exists $cache_file]} {
            # FAST PATH: Read directly from disk, bypass Python
            set fp [open $cache_file r]
            set ::cached_frustration_data [read $fp]
            close $fp
            
            # CRITICAL: Sync trackers so we don't immediately trigger a recalculation!
            set ::last_physics_pos $curr_phys
            set ::last_visual_pos $curr_vis
            set needs_redraw 1
        }
    }

    # -----------------------------------------------------------------
    # TIER 1: Did the physics change? (Cache miss OR interactive drag)
    # -----------------------------------------------------------------
    if {$curr_phys ne $::last_physics_pos} {
        
        set ::last_physics_pos $curr_phys
        set ::last_visual_pos $curr_vis
        
        # Send to Python
        puts $::py_pipe [join $curr_phys " "]
        flush $::py_pipe
        
        set new_data ""
        while {1} {
            set char_count [gets $::py_pipe line]
            if {$char_count == -1} {
                puts "ERROR: python script crashed!"
                catch {close $::py_pipe} 
                return
            } 
            if {$line eq "__END_OF_CALC__"} { break } 
            append new_data $line "\n"
        } 
        
        set ::cached_frustration_data $new_data
        set needs_redraw 1

        # WRITE CACHE: Save calculation to disk for instant loading next time.
        # If this was an interactive drag, it overwrites the stale frame cache.
        if {[catch {
            set fp [open $cache_file w]
            puts -nonewline $fp $new_data
            close $fp
        } err]} {
            puts "Warning: Could not write cache file $cache_file ($err)"
        }

    # -----------------------------------------------------------------
    # TIER 2: Did only visual anchors change? (e.g., smoothing window)
    # -----------------------------------------------------------------
    } elseif {$curr_vis ne $::last_visual_pos} {
        set ::last_visual_pos $curr_vis
        set needs_redraw 1
    }

    # Redraw only if something actually moved or frame changed
    if {$needs_redraw} {
        draw_lines_from_cache $molid
    }

    set ::coord_poll_id [after 10 auto_update_on_move]
}

# Establish connection to the persistent Python script
#set ::py_pipe [open "|python3 interactive_update.py" r+]
#fconfigure $::py_pipe -buffering line -blocking 1

# Global state trackers
set ::active_molid -1
set ::last_physics_pos ""
set ::last_visual_pos ""
set ::cached_frustration_data ""
array set ::mol_initial_atoms {}


# Setup persistent cache directory in user's home folder
set ::frust_cache_dir [file join ~ ".interactive-frustratometer" "frustratograms"]
file mkdir $::frust_cache_dir

# Add a tracker for the frame number
set ::last_frame -1

# main loop
catch {after cancel $::coord_poll_id}
auto_update_on_move

