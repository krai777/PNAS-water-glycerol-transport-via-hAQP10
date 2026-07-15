lassign $argv psf dcd
set mol [mol new ../../hAQP10_wi.psf type psf waitfor all]
mol addfile ../wrap_JAN21.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set numframe [molinfo top get numframes]
puts "Total frames: $numframe"
package require pbctools

# List of indices to analyze
set indices {161484 161806 161848 161344}

# Loop through each index
foreach idx $indices {
    puts "Processing index $idx"
    set out [open "wat_arnd_${idx}.dat" w]
    
    # Loop through each frame
    for {set i 0} {$i < $numframe} {incr i} {
        animate goto $i
        set wat_arnd [atomselect top "name OH2 and pbwithin 5 of index $idx"]
        set number [$wat_arnd num]
        puts $out "$i $number"
        $wat_arnd delete
    }
    
    close $out
    puts "Completed analysis for index $idx"
}

puts "Analysis complete for all indices"
quit
