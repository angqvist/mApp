res
set term postscript eps enh col "Helvetica" 22
set out "eps/fit.eps"
set enc iso_8859_1

#--------------------------------------------------------
col_red    = "#D7191C"
col_orange = "#FDAE61"
col_yellow = "#FFFFBF"
col_blue1  = "#ABD9E9"
col_blue2  = "#2C7BB6"
#------------
set style incr user
LW=6
LS=0
PS = 1*2.4
LS=LS+1 ; set style line LS lt 1 pt  4 ps PS lw LW lc rgb col_red
LS=LS+1 ; set style line LS lt 2 pt  6 ps PS lw LW lc rgb col_orange
LS=LS+1 ; set style line LS lt 3 pt  8 ps PS lw LW lc rgb col_blue2
LS=LS+1 ; set style line LS lt 3 pt  8 ps PS lw LW lc rgb col_blue2
LS=LS+1 ; set style line LS lt 1 pt  4 ps PS lw 10 lc rgb "#bbbbbb"
LS=LS+1 ; set style line LS lt 1 pt  6 ps PS lw 10 lc rgb "#bbbbbb"
LS=LS+1 ; set style line LS lt 1 pt  8 ps PS lw 10 lc rgb "#bbbbbb"
LS=LS+1 ; set style line LS lt 4 pt 10 ps PS lw LW lc rgb col_blue1
LS=10
LS=LS+1 ; set style line LS lt 1 lw 8 lc rgb "#bbbbbb"
LS=LS+1 ; set style line LS lt 2 lw 8 lc rgb "#bbbbbb"
#----
LW=4
LS=10
set style arr 11 lt 1 lw LW lc 0 filled head
set style arr 12 lt 1 lw LW lc 0 filled heads
set style arr 13 lt 3 lw LW lc 0 nohead

#--------------------------------------------------------
set xla "{/=24 Temperature}"
#set xra [0:1200]
#set xti 0,200
#----
set yla "{/=24 Energy (eV/unit cell) }"
#set yti 1e10
#set yra [10.83:11.5]
#set form y "%3.1e"
#----
set key Left rev spac 1.7 noautoti
set key left center at gr 0.05,0.35

x0 = 967 ; y0 = 0.18
#set arr from x0,gr y0 to x0,gr 0.01 as 11 front
#set lab "{/Helvetica-Oblique T_m}" at x0,gr y0+0.05 center front

set style fill solid 0.25

str = "%d{/Helvetica-Oblique %s}"

x1 = 900 ; x2 = 1100
d = 1 ; y1 = 10.84 ; y2 = 10.853
#set obj rect from x1,y1 to x2,y2 fc rgb col_red behind fs noborder solid 0.1
set obj rect from x1,y1 to x2,y2 fs empty border 1 lw 6 front



p \
  "ET_55p.data" u 1:2 w p,\
  "EV_55p.data" u 1:2 w p,\
  x w l