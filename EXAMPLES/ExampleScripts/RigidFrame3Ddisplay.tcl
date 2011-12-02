set displayType "PLAN"
#set displayType "PERSPECTIVE"
	
# a window showing the displaced shape
recorder display g3 10 10 300 300 -wipe

if {$displayType == "PERSPECTIVE"} {
	prp -7500 -5000 50000
	vrp 0 -500 250
	vup 0 0 1
	vpn 0 -1 0
	viewWindow -200 400 -300 300
}

if {$displayType == "PLAN"} {
	prp 0 0 1000
	vrp 0 0 0
	vup 0 -1 0
	vpn 0 0 -1
	viewWindow -200 200 -200 200
}

plane 0 1e12
port -1 1 -1 1
projection 1
fill 0
display 1 0 10