from packages import *

# LIST OF FUNCTIONS

def read_file_adv(file,index,r_1,r_2):
    data=np.loadtxt(file)
    return data[r_1:r_2,index]

def read_file(file,index):
    data=np.loadtxt(file)
    return data[:,index]

def linregress(x_ax,y_ax):
    m,c,r,p,s=sts.linregress(x_ax,y_ax)
    return m,s

def read_pert():
    path=file_dir+'energy_vs_time.dat'
    rtc[:]=read_file_adv(path,5,0,T)

def time_stamp():
    now    = datetime.now() # current date and time
    stamp  = now.strftime("%H%M%S_%m%d%Y")
    return stamp

# rc('text', usetex=True)
# rc('font',**{'family':'serif','serif':['Times New Roman']})
# rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans']})
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

# ax.xaxis.set_major_locator(FixedLocator(x_ticks))
# ax3.xaxis.set_major_locator(LogLocator(base=10))
# ax2.yaxis.set_label_coords(1.05,0.68)

# ax.text(0.82,0.055, r'$ \mathsf{Burgers}$', fontsize=4,color=plCol3)
# rect = patches.Rectangle((0.69, 0.042), 0.48, 0.017, linewidth=0.15, edgecolor=plCol3, facecolor='none')
# ax.add_patch(rect)

#  LINE STYLES
# linestyle_tuple = [
#      ('loosely dotted',        (0, (1, 10))),
#      ('dotted',                (0, (1, 1))),
#      ('densely dotted',        (0, (1, 1))),
#
#      ('loosely dashed',        (0, (5, 10))),
#      ('dashed',                (0, (5, 5))),
#      ('densely dashed',        (0, (5, 1))),
#
#      ('loosely dashdotted',    (0, (3, 10, 1, 10))),
#      ('dashdotted',            (0, (3, 5, 1, 5))),
#      ('densely dashdotted',    (0, (3, 1, 1, 1))),
#
#      ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
#      ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
#      ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]

# MARKERS
# "." 	m00 	point
# "," 	m01 	pixel
# "o" 	m02 	circle
# "v" 	m03 	triangle_down
# "^" 	m04 	triangle_up
# "<" 	m05 	triangle_left
# ">" 	m06 	triangle_right
# "1" 	m07 	tri_down
# "2" 	m08 	tri_up
# "3" 	m09 	tri_left
# "4" 	m10 	tri_right
# "8" 	m11 	octagon
# "s" 	m12 	square
# "p" 	m13 	pentagon
# "P" 	m23 	plus (filled)
# "*" 	m14 	star
# "h" 	m15 	hexagon1
# "H" 	m16 	hexagon2
# "+" 	m17 	plus
# "x" 	m18 	x
# "X" 	m24 	x (filled)
# "D" 	m19 	diamond
# "d" 	m20 	thin_diamond
# "|" 	m21 	vline
# "_" 	m22 	hline
