# Load PSF structure file
assign $argv psf dcd
mol new ../../hAQP10_wi.psf type psf
mol addfile ../wrap_JAN21.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set numframe [molinfo top get numframes]
puts "Total frames: $numframe"
package require pbctools

# Define the analysis parameters (modified to only include residue ID and chain ID)
set analysis_params {
    {161806 B}
    {161848 C}
    {161484 D}
    {161344 D}
    # Add more entries as needed in the format: {residue_id chain_id}
}

# Load molecule and process each analysis set
foreach params $analysis_params {
    # Extract parameters
    lassign $params residue_id chain_id
    
    # Create output filename
    set outfile [open "glyangle_${residue_id}_${chain_id}.dat" w]
    puts $outfile "Frame  Z_Position  C1_C3_Angle"
    
    # Select the protein (pore structure)
    set prot [atomselect top "chain $chain_id"]
    
    # Process all frames
    for {set i 0} {$i < $numframe} {incr i} {
        animate goto $i
        
        # Update selections
        $prot update
        # Define pore axis (Z-direction)
        set pore_axis {0 0 1}
        
        # Select glycerol atoms (C1 and C3)
        set c1 [atomselect top "(same residue as index $residue_id) and name C1"]
        set c3 [atomselect top "(same residue as index $residue_id) and name C3"]
        
        # Ensure selection exists
        if { [$c1 num] == 0 || [$c3 num] == 0 } {
            continue
        }
        
        # Get coordinates
        set pos_c1 [lindex [$c1 get {x y z}] 0]
        set pos_c3 [lindex [$c3 get {x y z}] 0]
        
        # Compute vector from C1 to C3
        set vec_C1_C3 [vecsub $pos_c3 $pos_c1]
        
        # Compute dot product between C1→C3 vector and pore axis
        set dot [vecdot $vec_C1_C3 $pore_axis]
        
        # Compute magnitudes
        set mag_C1_C3 [veclength $vec_C1_C3]
        set mag_pore_axis [veclength $pore_axis]
        
        # Compute angle (in degrees)
        set angle [expr acos($dot / ($mag_C1_C3 * $mag_pore_axis)) * 180 / 3.14159265]
        
        # Compute glycerol center for Z-position
        set gly_center [measure center [atomselect top "same residue as index $residue_id"] weight mass]
        set gly_z [lindex $gly_center 2]
        
        # Write data to file
        puts $outfile "$i  $gly_z  $angle"
        
        # Free selections
        $c1 delete
        $c3 delete
    }
    
    close $outfile
    puts "Completed analysis for residue $residue_id with chain $chain_id"
}

puts "All calculations completed."
quit
