set title "CapDEM"
set xlabel "Saturation"
set ylabel "Water Potential"
set xrange [0:1]
set yrange [-5:30]
set grid

plot "data/g1e-3_s0109/ahistory" every 2::2 using 7:8 title "g1e-3_s0109" w p,\
"data/g1e-3_s0208/ahistory" every 2::2 using 7:8 title "g1e-3_s0208" w p,\
"data/g1e-3_s0406/ahistory" every 2::2 using 7:8 title "g1e-3_s0406" w p,\
"data/g1e-5_s0109/ahistory" every 2::2 using 7:8 title "g1e-5_s0109" w l,\
"data/g1e-5_s0208/ahistory" every 2::2 using 7:8 title "g1e-5_s0208" w l,\
"data/g1e-5_s0406/ahistory" every 2::2 using 7:8 title "g1e-5_s0406" w l,\
"data/g1e-8_s0109/ahistory" every 2::2 using 7:8 title "g1e-8_s0109" w d,\
"data/g1e-8_s0208/ahistory" every 2::2 using 7:8 title "g1e-8_s0208" w d,\
"data/g1e-8_s0406/ahistory" every 2::2 using 7:8 title "g1e-8_s0406" w d
