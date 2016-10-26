sed -i "s/_create(\(.*\)8UL/_create(\1$1UL/" nbody.Pthreads.c 
#s/LLmrf_write_begin(PNchannel_\([a-zA-Z-_]*\),/printf(\"wb{\1}\\n\"); LLmrf_write_begin(PNchannel_\1,/g
#s/LLmrf_write_end(PNchannel_\([a-zA-Z-_]*\),/printf(\"we{\1}\\n\"); LLmrf_write_end(PNchannel_\1,/g 
