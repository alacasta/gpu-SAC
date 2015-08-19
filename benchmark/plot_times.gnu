    set terminal pngcairo dashed size 900.0, 300.0 enhanced font 'Corbel,10'
    set out "cpuVSgpu.png"
    set key below
	unset key
	set xtics nomirror tc rgb "#808080"
	set ytics nomirror tc rgb "#808080"	
    set ylabel "Time (s)" tc rgb "#808080"	
    set xlabel "Number of elements" tc rgb "#808080"	
	set style line 103 lc rgb'#808080' lt 2 lw 2
    set style line 102 lc rgb'#808080' lt 0 lw 1
	set logscale y
	set logscale x	
	#set y2range[0.0:2.0]
    load '../palettes/diverging/Spectral.plt'	
	set palette negative
	set format x "10^{%L}"
	set format y "10^{%L}"
	set mxtics 10
	set grid mxtics mytics xtics ytics lt 1 lc rgb 'gray70', lt 1 lc rgb 'gray90'
	set key below
    plot "sac_cpu_times.out" u ($1*$1):2 w lp lt 1 lw 2 lc  rgb'#66ccff' t'1-Core CPU',\
	"sac_gpu_times.out" u ($1*$1):2 w lp lt 1 lw 2 lc  rgb'#00CC66' t 'GPU'