To Compile - goto the water_nsquared_transactional/ folder and type : 
make clean; make water 
an executable , hash, would be created in the same folder

The number of threads could be changed, once you browse through the code 
and you could change the granularity and play with the barriers.

It scales well(better than the locks version)  from 1 to 2 to 4 threads and then does not after 4. 
The execution times decrease till 4 threads then increase as you can see in the file performance.xls 
<look at the percentages also>
