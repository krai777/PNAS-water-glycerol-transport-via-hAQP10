lassign $argv psf dcd
mol load psf $psf dcd $dcd
set numframe [molinfo top get numframes]
package require pbctools
set out [open "watnum.dat" w ]      
for {set i 0} {$i <$numframe} {incr i} {
    animate goto $i 
    set z1 [molinfo 0 get c]
    set cutoff [expr ($z1) +15]
    set o1 [atomselect top "(abs(z) > $cutoff and name OH2)"]
    set s1 [atomselect top "(abs(z) > $cutoff and name SOD)"]
    set c1 [atomselect top "(abs(z) > $cutoff and name CLA)"]
    set g1 [atomselect top "(abs(z) > $cutoff and resname MGLY and name C2)"]
    set o2 [atomselect top "(z<-15 and segname WATUP) or (segname WATDOWN and z>15)"]
    set s2 [atomselect top "(z<-15 and segname SODUP) or (segname SODDOWN and z>15)"]
    set c2 [atomselect top "(z<-15 and segname CLAUP) or (segname CLADOWN and z>15)"]
    set g2 [atomselect top "(z<-15 and segname GLYUP) or (segname GLYDOWN and z>15)"]
    set z_list [$o1 get z]
    set no1 0 
    foreach w $z_list {
        set l [expr abs(int($w)/int($cutoff))]
        set no1 [expr $l+$no1]
    }  
    set z_list_s [$s1 get z]
    set ns1 0
    foreach w $z_list_s {
        set l [expr abs(int($w)/int($cutoff))]
        set ns1 [expr $l+$ns1]
    }
    set z_list_c [$c1 get z]
    set nc1 0
    foreach w $z_list_c {
        set l [expr abs(int($w)/int($cutoff))]
        set nc1 [expr $l+$nc1]
    }
    set z_list_g [$g1 get z]
    set ng1 0
    foreach w $z_list_g {
        set l [expr abs(int($w)/int($cutoff))]
        set ng1 [expr $l+$ng1]
    }
    set no2 [$o2 num]
    set ns2 [$s2 num]
    set nc2 [$c2 num]
    set ng2 [$g2 num]
    set no [expr $no2 +$no1]
    set ns [expr $ns2 +$ns1]
    set nc [expr $nc2 +$nc1]
    set ng [expr $ng2 +$ng1]
    puts $out "[expr ($i+1) ] $no $ns $nc $ng"
}
close $out
quit
