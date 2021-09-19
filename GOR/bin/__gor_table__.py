#!/usr/bin/python
#!python


# GOR_Tables

# GOR Table_1 [ Alpha Helix Conformation ] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
helix_table = {"G":[-5,-10,-15,-20,-30,-40,-50,-60,-86,-60,-50,-40,-30,-20,-15,-10,-5],
	 "A":[5,10,15,20,30,40,50,60,65,60,50,40,30,20,15,10,5],
	 "V":[0,0,0,0,0,0,5,10,14,10,5,0,0,0,0,0,0],
	 "L":[0,5,10,15,20,25,28,30,32,30,28,25,20,15,10,5,0],
	 "I":[5,10,15,20,25,20,15,10,6,0,-10,-15,-20,-25,-20,-10,-5],
	 "S":[0,-5,-10,-15,-20,-25,-30,-35,-39,-35,-30,-25,-20,-15,-10,-5,0],
	 "T":[0,0,0,-5,-10,-15,-20,-25,-26,-25,-20,-15,-10,-5,0,0,0],
	 "D":[0,-5,-10,-15,-20,-15,-10,0,5,10,15,20,20,20,15,10,5],
	 "E":[0,0,0,0,10,20,60,70,78,78,78,78,78,70,60,40,20],
	 "N":[0,0,0,0,-10,-20,-30,-40,-51,-40,-30,-20,-10,0,0,0,0],
	 "Q":[0,0,0,0,5,10,20,20,10,-10,-20,-20,-10,-5,0,0,0],
	 "K":[20,40,50,55,60,60,50,30,23,10,5,0,0,0,0,0,0],
	 "H":[10,20,30,40,50,50,50,30,12,-20,-10,0,0,0,0,0,0],
	 "R":[0,0,0,0,0,0,0,0,-9,-15,-20,-30,-40,-50,-50,-30,-10],
	 "F":[0,0,0,0,0,5,10,15,16,15,10,5,0,0,0,0,0],
	 "Y":[-5,-10,-15,-20,-25,-30,-35,-40,-45,-40,-35,-30,-25,-20,-15,-10,-5],
	 "W":[-10,-20,-40,-50,-50,-10,0,10,12,10,0,-10,-50,-50,-40,-20,-10],
	 "C":[0,0,0,0,0,0,-5,-10,-13,-10,-5,0,0,0,0,0,0],
	 "M":[10,20,25,30,35,40,45,50,53,50,45,40,35,30,25,20,10],
	 "P":[-10,-20,-40,-60,-80,-100,-120,-140,-77,-60,-30,-20,-10,0,0,0,0] }

# GOR Table_2 [ Extended Conformation ] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
extended_table = {"G":[10,20,30,40,40,20,0,-20,-42,-20,0,20,40,40,30,20,-10],
	 "A":[0,0,0,0,-5,-10,-15,-20,-23,-20,-15,-10,-5,0,0,0,0],
	 "V":[0,0,-10,-20,0,20,40,60,68,60,40,20,0,-20,-10,0,0],
	 "L":[0,0,0,0,0,5,10,20,23,20,10,5,0,0,0,0,0],
	 "I":[0,-10,-20,-10,0,20,40,60,67,60,40,20,0,-10,-20,-10,0],
	 "S":[0,10,20,10,0,-5,-10,-15,-17,-15,-10,-5,0,10,20,10,0],
	 "T":[5,10,15,20,15,15,10,10,13,10,10,15,15,20,15,10,5],
	 "D":[0,5,10,15,20,0,-20,-30,-44,-30,-20,0,0,0,0,0,0],
	 "E":[-10,-15,-20,-25,-30,-35,-40,-45,-50,-55,-60,60,-50,-40,-30,-20,-10],
	 "N":[10,30,50,30,20,0,-15,-30,-41,-30,-15,0,20,30,50,30,10],
	 "Q":[0,0,0,0,0,-5,-10,0,12,20,30,40,50,50,40,30,15],
	 "K":[-5,-10,-15,-20,-30,-40,-50,-40,-33,-20,-10,0,10,10,0,0,0],
	 "H":[-10,-20,-40,-20,-10,0,-10,-20,-25,-35,-30,-25,-20,-15,-10,-5,0],
	 "R":[0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0],
	 "F":[0,0,0,0,0,5,10,20,26,10,-10,-30,-60,-65,-60,-40,-20],
	 "Y":[0,5,10,15,20,25,30,35,40,35,30,25,20,15,10,5,0],
	 "W":[0,0,0,0,0,-10,-10,-10,-10,-10,-10,-15,-20,-25,-30,-20,-10],
	 "C":[0,0,0,0,0,10,20,30,44,30,20,10,0,0,0,0,0],
	 "M":[-10,-20,-30,-40,-40,-30,0,10,23,10,0,-30,-40,-40,-30,-20,-10],
	 "P":[10,20,30,30,20,10,0,-10,-18,-20,-10,10,30,40,30,20,10] }

# GOR Table_3 [ Turns Conformation ] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
turns_table = {"G":[0,0,0,0,10,30,55,55,57,40,0,0,0,0,0,0,0],
	 "A":[0,0,0,-10,-20,-30,-40,-50,-50,-40,-30,-20,-10,0,0,0,0],
	 "V":[0,0,0,0,-10,-20,-30,-40,-60,-40,-30,-20,-10,0,0,0,0],
	 "L":[0,0,0,-10,-20,-30,-40,-50,-56,-20,-10,0,0,0,0,0,0],
	 "I":[0,0,0,0,0,-10,-20,-30,-46,-40,-10,0,0,20,30,20,10],
	 "S":[0,-10,-20,-20,10,15,20,25,26,25,20,15,10,0,0,0,0],
	 "T":[0,10,20,20,20,15,18,5,3,5,10,15,20,20,20,10,0],
	 "D":[0,0,0,0,0,0,5,10,31,10,5,0,0,0,0,0,0],
	 "E":[0,-5,-10,-15,-20,-30,-40,-45,-47,-20,0,10,5,0,0,0,0],
	 "N":[0,0,0,10,20,30,35,40,42,40,35,30,20,10,5,0,0],
	 "Q":[10,20,30,25,20,15,10,5,4,20,30,40,50,60,50,40,20],
	 "K":[-10,-20,-30,-40,-25,-10,0,10,10,10,0,-20,-30,-20,-10,-5,0],
	 "H":[0,0,0,0,0,0,0,0,-3,0,10,20,30,20,10,0,0],
	 "R":[0,0,0,0,0,0,0,10,21,30,40,30,20,10,0,0,0],
	 "F":[0,0,0,0,0,-5,-10,-15,-18,-15,0,15,30,25,20,10,0],
	 "Y":[0,0,0,5,15,15,20,25,29,25,20,15,15,5,0,0,0],
	 "W":[0,0,0,10,20,30,40,80,36,-30,30,40,50,60,70,40,20],
	 "C":[20,40,50,60,60,55,50,45,44,40,35,30,25,20,15,10,5],
	 "M":[-5,-15,-20,-25,-30,-35,-40,-45,-48,-45,-40,-35,-30,-25,-20,-15,-5],
	 "P":[10,20,30,40,50,70,10,-90,36,90,10,0,0,0,0,0,0] }

# GOR Table_4 [ Coil Conformation ] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
coil_table = {"G":[0,0,0,0,10,30,40,45,49,45,40,30,10,0,0,0,0],
	 "A":[0,0,0,0,-5,-10,-20,-25,-25,-25,-20,-15,-10,-5,0,0,0],
	 "V":[0,0,0,0,-10,-20,-25,-30,-35,-30,-25,-20,-10,0,0,0,0],
	 "L":[0,0,0,-10,-20,-30,-40,-30,-20,-20,-10,0,0,0,0,0,0],
	 "I":[0,0,0,0,0,-10,-20,-30,-33,-30,-10,0,10,20,30,20,0],
	 "S":[0,-10,-20,-20,10,15,20,25,50,25,20,15,10,0,0,0,0],
	 "T":[0,10,20,30,20,15,10,15,17,15,10,15,20,30,20,10,0],
	 "D":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
	 "E":[0,0,10,20,40,20,0,-10,-44,-40,-20,-10,0,0,0,0,0],
	 "N":[0,0,0,10,20,30,35,40,46,40,35,30,20,10,0,0,0],
	 "Q":[10,20,30,25,20,15,10,0,-5,20,30,40,50,60,50,40,20],
	 "K":[-10,-20,-30,-40,-25,-20,-10,-8,-8,0,0,-20,-30,-20,-10,-5,0],
	 "H":[0,0,0,0,0,0,0,10,16,15,10,10,10,10,5,0,0],
	 "R":[0,0,0,0,0,0,0,0,-12,0,20,30,20,10,0,0,0],
	 "F":[0,0,0,0,0,-5,-10,-20,-41,-20,0,15,30,25,20,10,0],
	 "Y":[0,0,0,0,0,0,0,0,-6,0,0,0,0,0,0,0,0],
	 "W":[0,0,0,10,20,30,40,20,12,20,30,40,50,60,70,40,20],
	 "C":[0,0,0,0,0,0,-10,-30,-47,-30,-10,0,0,0,0,0,0],
	 "M":[0,-5,-10,-15,-20,-25,-30,-40,-41,-40,-30,-25,-20,-15,-10,-5,0],
	 "P":[0,0,10,20,30,40,50,55,58,50,10,0,0,0,0,0,0] }
