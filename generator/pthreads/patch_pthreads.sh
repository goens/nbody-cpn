if [ "$#" -ne 2 ]; then
  echo "Illegal number of parameters"
  echo "Usage: ./patch_pthreads.h. <target_file> <new_channel_size>"
  exit 1 
fi

sed -i "s/_create(\(.*\)8UL/_create(\1$2UL/" $1
#s/LLmrf_write_begin(PNchannel_\([a-zA-Z-_]*\),/printf(\"wb{\1}\\n\"); LLmrf_write_begin(PNchannel_\1,/g
#s/LLmrf_write_end(PNchannel_\([a-zA-Z-_]*\),/printf(\"we{\1}\\n\"); LLmrf_write_end(PNchannel_\1,/g 
