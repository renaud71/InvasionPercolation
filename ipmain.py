
from ip import IP	
import time

start_time = time.time()
angles = [0,15,30,45,60] #vinkler


for i in range(len(angles)):

	c=IP(Ncol=344, Nrow=260, theta=angles[i])
	c.invasion()
	print("--- %s seconds ---" % (time.time() - start_time))