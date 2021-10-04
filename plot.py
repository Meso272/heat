import os

start=20000
end=25001
step=100
for i in range(start,end,step):
	path="mpires/%d.dat" % i
	command="PlotSliceImage -f -i %s -2 200 200 -m INDV -n ORI -o %s.png" % (path,path)
	os.system(command)