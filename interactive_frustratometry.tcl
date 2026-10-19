# ABOUT THIS SCRIPT 

# USE WITH VMD TO CALCULATE FRUSTRATION ON-THE-FLY 
# AS COORDINATES OF A STRUCTURE ARE CHANGED.
#
# IF MULTIPLE MOLECULES ARE LOADED,
# THE FRUSTRATOGRAM WILL BE CALCULATED FOR THE TOP MOLECULE
# AND RE-CALCULATED AS THE TOP MOLECULE IS CHANGED.
# THIS REQUIRES CLOSING AND RE-LAUNCHING THE PYTHON PIPE,
# WHICH MAY CAUSE A DELAY OF SEVERAL SECONDS.
#
# ONCE THE PYTHON PIPE IS INITIALIZED, RECALCULATION OF THE
# FRUSTRATOGRAM FOLLOWING AN ADJUSTMENT OF THE MOLECULAR COORDINATES
# IS OFTEN FAST ENOUGH TO KEEP UP WITH A MOVIE
# PLAYED AT THE NORMAL FRAME RATE.
#
# YOU MAY BE PROMPTED FOR THE SEQUENCE AND LOCATIONS OF THE
# CHAIN BREAKS IN YOUR STRUCTURE IF THE SEQUENCE CANNOT BE
# DETERMINED FROM YOUR STRUCTURE FILE (E.G., IN THE CASE
# THAT YOUR STRUCTURE FILE USES IGL/IPR/NGP RESIDUE NAMES
# INSTEAD OF STANDARD 3-LETTER RESIDUE NAMES)

# =====================================================================

# PROCEDURE DEFINITIONS

# procedure to clear the lines of the frustratogram, 
# as stored in the global frustratogram_lines array
proc clear_frustratogram_lines {molid} {
    foreach gid $::frustratogram_lines {
        catch {graphics $molid delete $gid}
    }
    set ::frustratogram_lines {}
}

# Tkinter Pop-up for User Input
# (need to ask the user for the sequence and chain breaks if
#  first attempt at structure initialization fails)
proc ask_user_for_python {prompt_text} {
    package require Tk
    set ::user_response ""

    set t [toplevel .py_prompt]
    wm title $t "Frustratometer: Input Required"
    wm attributes $t -topmost 1 ;# put the window on top of VMD

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
# procedure to terminate this script upon failure of the python 
# initialization code found in start_python_for_mol
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
    clear_frustratogram_lines $::active_molid
}


# =====================================================================
# procedure to launch the python script, which will 
# initialize the frustration-calculating python object (such as AWSEM)
# from the appropriate structure file on disk
# =====================================================================
proc start_python_for_mol {molid} {
    # close the existing python pipe (if it exists)
    if {[info exists ::py_pipe]} { catch {close $::py_pipe} }

    # Extract the original file path.
    # 'molinfo get filename' returns a nested list: {{file1.pdb} {file2.dcd}}
    # We grab the very first file, which should be a PDB that is compatible
    # with the structure parser in the main python frustratometer code
    set struct_file [lindex [lindex [molinfo $molid get filename] 0] 0]

    # Verify that the file exists in the location that we expect
    if {![file exists $struct_file]} {
        kill_frustration_script "Cannot locate original structure file on disk: $struct_file"
        return 0
    }

    puts "Starting Python server for molecule $molid using file: $struct_file"
    # launch the update_frustration.py script and set up pipe to its stdin/stdout
    if {[catch {open "|python3 -u update_frustration.py \"$struct_file\"" r+} ::py_pipe]} {
        kill_frustration_script "pipe failure: $::py_pipe"
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
            kill_frustration_script "Python failed to start.\nPython error was:\n$python_error"
            return 0
        }

        if {$line eq "__READY__"} {
            puts "Frustratometer initialized and ready!"
            break
        } elseif {$line eq "__PROMPT__"} {
            # Python wants input. 
            # The very next line in the pipe is the question for the user
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

    # initialize or clear preexisting coorindates and data 
    # that are no longer relevant
    set ::last_physics_pos "PLACEHOLDER"
    set ::last_visual_pos "PLACEHOLDER"
    set ::cached_frustration_data "PLACEHOLDER"

    return 1
}


# =====================================================================
#    procedure to draw the frustratogram lines
#
#    Whenever a CA atom is moved, we need to at least update the 
#    endpoints of some of the frustratogram lines on the screen.
#    
#    We will likely also need to send the new coordinates to python 
#    and receive new frustration data (which chains/residues to 
#    connect with lines and what kinds of interactions those lines 
#    describe (short, long, water=mediated; minimally frustrated,
#    highly frustrated)) back from python. This is handled by 
#    the auto_update_on_move procedure.
#
#    Here, we just clear stale data and then draw lines based on the 
#    new frustration data from the global scope. 
#
# =====================================================================
proc draw_lines_from_cache {molid} {

    # clear any previously drawn lines
    clear_frustratogram_lines $molid

    # In rare cases, python might not identify any minimally or highly 
    # frustrated contacts. Since we've already cleared any preexisting 
    # lines, the correct thing to do is to return. But this 
    # circumstance was probably not expected by the user, so we should 
    # explain what is going on.
    if {$::cached_frustration_data eq ""} { 
        puts "Python found no minimally or highly frustrated contacts, so there is nothing to draw!"    
        return
    }

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

        # determine whether/how to draw line;
        # "neutral" pairs are generally not expected to be 
        # returned through the pipe
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

        # draw the line and append its id to our list of graphics 
        # so that it can be deleted later upon a change of coordinates 
        # or switching to a new top molecule
        if {[info exists coords($c1,$resid1)] && [info exists coords($c2,$resid2)]} {
            graphics $molid color $color
            set gid [graphics $molid line $coords($c1,$resid1) $coords($c2,$resid2) style $style width 2]
            lappend ::frustratogram_lines $gid
        }
    }
}

# =====================================================================
#  procedure that periodically checks whether coordinates have changed 
#  and moves the locations of the ends of the lines drawn on the screen,
#  first recalculating which lines need to be drawn, if necessary.
# =====================================================================
proc auto_update_on_move {} {

    # select the index of the top molecule 
    set top_mol [molinfo top]

    # if the top molecule changed (including the case that we deleted
    # the only open molecule, in which VMD sets top to -1),
    # then delete the frustratogram lines 
    if {$top_mol != $::active_molid} {
        clear_frustratogram_lines $::active_molid
    }

    # A top index of -1 means that there are no molecules open in VMD.
    # In this case, there is no calculation to run, so we 
    # schedule another run of this procedure a few ms in the future
    # and then exit.
    if {$top_mol == -1} {
        set ::coord_poll_id [after 10 auto_update_on_move]
        return
    }

    # check if the user changed the top molecule 
    # and, if so, reboot the frustration-calculating python code
    # by calling start_python_for_mol
    if {$top_mol != $::active_molid} {
        if {![start_python_for_mol $top_mol]} {
            # If the python startup failed, then 
            # return now and don't schedule another run 
            # of this procedure. The user will need to 
            # re-source this script with a different
            # molecule selected as top if they want to 
            # continue to display frustratograms.
            return
        }
    }

    # set the active molecule to the top molecule
    set ::active_molid $top_mol

    # get the coordinates of the CB atoms (CA in the case of glycine), 
    # which determine the pairwise interaction energies
    set sel_phys [atomselect $::active_molid "name CB or (resname GLY and name CA)"]
    set curr_phys [$sel_phys get {x y z}]
    $sel_phys delete

    # get the coordinates of the CA atoms, which will become the  
    # endpoints of the frustratogram lines that we draw
    set sel_vis [atomselect $::active_molid "name CA"]
    set curr_vis [$sel_vis get {x y z}]
    $sel_vis delete

    # If either of these arrays (which are just strings of different lengths in TCL)
    # changed since last time ("last time" being an empty string if this is the 
    # first time through the loop since we most recently set/changed the top molecule),
    # then we will recalculate the frustratogram (if needed) 
    # and update its visual representation.
    set needs_redraw 0
    if {$curr_phys ne $::last_physics_pos} {

        set ::last_physics_pos $curr_phys
        set ::last_visual_pos $curr_vis

        # Flatten the coordinate list { {x1 y1 z1} {x2 y2 z2} } -> "x1 y1 z1 x2 y2 z2"
        # and send it straight down the pipe to the python code
        puts $::py_pipe [join $curr_phys " "]
        flush $::py_pipe

        # Read instructions (in the form of a long string)
        # from the pipe until we reach the special 
        # __END_of_CALC__ signal
        set new_data ""
        while {1} {
            set char_count [gets $::py_pipe line]

            # check if our pipe broke due to the python script crashing
            if {$char_count == -1} {
                puts "ERROR: python frustration calculation crashed. Did you change the topology of your molecule?"
                # close the broken pipe
                catch {close $::py_pipe}
                # If the python frustration calculation failed, then 
                # return now and don't schedule another run 
                # of this procedure. It is likely that the user has 
                # edited the topology of the active molecule (not supported)
                # or there is an error with this script.
                return
            }

            if {$line eq "__END_OF_CALC__"} { break }
            append new_data $line "\n"
        }

        # update the data to be passed into the drawing function
        # draw_lines_from_cache a few lines below this point
        set ::cached_frustration_data $new_data
        set needs_redraw 1

    } elseif {$curr_vis ne $::last_visual_pos} {
        set ::last_visual_pos $curr_vis
        set needs_redraw 1
    }

    if {$needs_redraw} {
        draw_lines_from_cache $::active_molid
    }

    # Schedule this procedure to be called again in 10 ms 
    # so that the frustratograms updates quickly in response to 
    # changes in the coordinates resulting from changes in the 
    # trajectory frame or editing of the structure in the GUI.
    set ::coord_poll_id [after 10 auto_update_on_move]
}


# ======================================================================
# MAIN BLOCK


# active_molid is the id of the molecule that we just drew the 
# frustratogram for. Whenever active_molid differs from VMD's "top" 
# molecule, we update active_molid to the top molecule's id and 
# reboot the frustration-calculating python script, feeding it the 
# topology and coordinates of the new top molecule. When no  
# molecules are present in the VMD workspace, the active_molid is set to 
# -1, just like the value of "molinfo top". We set active_molid to -1
# when initializing this script as a placeholder. It will be checked 
# and updated immediately at the beginning of the auto_update_on_move
# procedure (called below).
set ::active_molid -1

# last_physics_pos holds the CB (CA in the case of glycine)
# coordinates of the molecular configuration used for the 
# most recent frustratogram draw. We check this periodically
# (once during every call to the auto_update_on_move procedure)
# to determine whether we need to recalculate the frustraogram
set ::last_physics_pos "PLACEHOLDER"

# last_visual_pos holds the CA coordinates that are to be 
# connected by frustratogram lines. We need to update the
# frustratogram by moving the lines (and very likely also
# recalculating the frustration) whenever these coordinates
# change.
set ::last_visual_pos "PLACEHOLDER"

# cached_frustration_data is a big string holding 
# concatenated lines of frustration information in the form 
# of frustratometer server output, meaning lines formatted like 
#     f'{resid1} {resid2} {chain1} {chain2} foo foo foo foo foo foo foo foo {style} {descriptor}\n'
# (the fields labeled "foo" don't matter for our purposes here;
# "style" is either "short", "long", or "water-mediated";
# and "descriptor" is either "minimally" or "highly")
set ::cached_frustration_data "PLACEHOLDER"

# frustratogram_lines holds the ids of the graphics (lines)
# associated with the frustratogram; 
# allows us to clean up the frustratogram without 
# destroying any other graphics the user might have created
set ::frustratogram_lines {}

# When this script is first launched, ::coord_poll_id is not defined,
# so $::coord_poll_id throws an error that is ignored because of the 
# catch command on the below line of code. However, if this script is
# sourced again, the main loop from the previous launch of this
# script is not automatically cancelled. The main loop, which is the
# body of the auto_update_on_move procedure, perpetuates itself by
# scheduling another invocation of itself for 10 ms in the future 
# using 'after 10 auto_update_on_move'. We can break the loop by 
# cancelling the next scheduled run of auto_update_on_move,
# whose command id is stored in coord_poll_id.
catch {after cancel $::coord_poll_id}

# The auto_update_on_move procedure implements the main loop. 
# If has the following parts:
#  - check which molecule we need to draw the frustratogram for 
#    (the top molecule)
#  - check whether the coordinates have been changed since the 
#    previous iteration of the loop
#  - recaculate the frustratogram (if needed)
#  - adjust the number color, and positions of the endpoints 
#    of the lines (if needed)
# Errors are expected to occur only within the python script 
# that parses the structures and calculates the frustration.
# When an error occurs in the python script, this TCL script 
# should print the error to the terminal that is running VMD 
# and then effectively terminate (choose not to schedule 
# another call to the auto_update_on_move function).
# Runtime errors that originate in this TCL script should be 
# viewed as bugs.
auto_update_on_move
