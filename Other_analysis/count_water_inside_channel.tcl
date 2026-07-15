mol load psf ../../hAQP10_wi.psf
mol addfile ../wrap_JAN21.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

set numframe [molinfo top get numframes]
package require pbctools

# Define chains
set chains {A B C D}

# Loop over chains
foreach chain $chains {
    # Open separate output file for each chain
    set outfile [open "wat_in_${chain}.dat" w]  ;# Output file for each chain

    for {set i 0} {$i < $numframe} {incr i} {
        animate goto $i  ;# Go to frame i

        # Measure the center of mass (COM) of the protein for the specific chain
        set prot [atomselect top "chain $chain"]
        set com_prot [measure center $prot weight mass]

        set com_x [lindex $com_prot 0]
        set com_y [lindex $com_prot 1]
        set com_z [lindex $com_prot 2]

        # Select water molecules based on the specific chain's COM
        set water_sel [atomselect top "water and noh and (x-$com_x)^2 + (y-$com_y)^2 < 81 and abs(z)<18"]

        set water_count [$water_sel num]
        puts "processed frame $i for chain $chain"

        # Write the frame number and water count to the corresponding file
        puts $outfile "$i\t$water_count"
    }

    close $outfile  ;# Close the output file for the current chain
}

quit
